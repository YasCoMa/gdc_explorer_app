import os
import sys
import json
import gzip
import subprocess
import shutil
import requests
import pandas as pd
import numpy as np

from tqdm import tqdm
import xml.etree.ElementTree as ET

from matplotlib import pyplot as plt


from sksurv.nonparametric import kaplan_meier_estimator

class DataWrangler:
    def __init__(self, fout = '../data_processed', cutoff_cases_expression_rnaseq = 2):
        self.cutoff_cases_expression_rnaseq = cutoff_cases_expression_rnaseq

        self.formats_datCategory = { "biospecimen": ["svs", "jpeg 2000"], "clinical": ["bcr xml"], "copy number variation": ["tsv", "txt"], "dna methylation": ["txt"], "proteome profiling": ["tsv"], "simple nucleotide variation": ["maf"], "transcriptome profiling": ["tsv"] }

        self.out = fout
        if( not os.path.exists(self.out) ):
            os.makedirs( fout )

    def _reduce_methylation_cg_mapping(self, filename, col_gene, col_chrom, platform, alldc):
        col_gene = int(col_gene)
        col_chrom = int(col_chrom)

        ds = {}
        f = open(filename)
        for line in f:
            if( line.startswith('cg') ):
                l = line.replace('\n','').split(',')

                cgid = l[0].split('_')[0]
                gene = l[col_gene]
                chrom = l[col_chrom].replace('chr', '')

                if( (gene != '') and (chrom != '') ):
                    ds[cgid] = { 'gene': gene, 'chromosome': chrom }
                    alldc[cgid] = { 'gene': gene, 'chromosome': chrom }
        f.close()

        opath = os.path.join( self.out, 'mapp_'+platform+'.json')
        json.dump( ds, open(opath, 'w') )

        return alldc

    def run_mapping_methylation(self):
        '''
        Calls:
        python3 process_data.py mapping_methylation /home/yasmmin/Downloads/in_met_mappings/illumina_humanmethylation27.csv 10 2 illumina_27

        python3 process_data.py mapping_methylation /home/yasmmin/Downloads/in_met_mappings/humanmethylation450_15017482_v1-2.csv 21 11 illumina_450

        python3 process_data.py mapping_methylation /home/yasmmin/Downloads/in_met_mappings/infinium-methylationepic-v-1-0-b5-manifest-file.csv 15 11 illumina_epic

        python3 process_data.py mapping_methylation /home/yasmmin/Downloads/in_met_mappings/EPIC-8v2-0_A1.csv 24 15 illumina_epicv2
        '''

        prev = 0
        opath = os.path.join( self.out, 'all_mapp.json')
        try:
            alldc = json.load( open(opath, 'r') )
            prev = len(alldc)
        except:
            alldc = {}

        filename = sys.argv[2]
        col_gene = sys.argv[3]
        col_chrom = sys.argv[4]
        platform = sys.argv[5]

        gdcnames_plat = ["illumina human methylation 27", "illumina human methylation 450", "illumina methylation epic", "illumina methylation epic v2"]
        plat_options = ["illumina_27", "illumina_450", "illumina_epic", "illumina_epicv2"]
        dc = dict( zip(plat_options, gdcnames_plat) )

        if( platform in plat_options ):
            alldc = self._reduce_methylation_cg_mapping(filename, col_gene, col_chrom, dc[platform].replace(' ','-'), alldc )

            opath = os.path.join( self.out, 'all_mapp.json')
            json.dump( alldc, open(opath, 'w') )
            if( len(alldc) > prev ):
                print(platform, ' - ', len(alldc))
        else:
            print("Please choose a valid platform for the methylation beads: illumina_27, illumina_450, illumina_epic or illumina_epicv2")

    def get_cases_info_by_project(self, project):
        basename = "%s_%s" %(project, 'cases')
        odir = os.path.join(self.out, "%s" %(basename) )
        if( not os.path.exists(odir) ):
            os.makedirs(odir)

        opath = os.path.join(odir, "cases_metadata.json")
        if( not os.path.isfile(opath) ):
            
            filters = { "op": "and",
                    "content":[
                        {
                            "op":"in",
                            "content":{
                                "field":"project.project_id",
                                    "value": [ project ]
                            }
                        }
                    ]
                }

            fields = [
                "case_id",
                "submitter_id",
                "samples.tumor_descriptor",
                "samples.tissue_type",
                "demographic.ethnicity",
                "demographic.gender",
                "demographic.race"
            ]

            fields = ",".join(fields)
            endpoint = "https://api.gdc.cancer.gov/cases"

            params = {
                "filters": json.dumps(filters),
                "fields": fields,
                "format": "json",
                "size": "2000"
            }
            response = requests.get(endpoint, params = params)
            d = response.json()
            
            df = {}
            for case in d['data']['hits']:
                if( 'demographic' in case ):
                    _id = case['case_id']
                    meta = case['demographic']
                    meta['submitter_id'] = case['submitter_id']
                    df[_id] = meta
            
            with open( opath, 'w') as f:
                json.dump( df, f)
        else:
            df = json.load( open(opath, 'r') )

        return df

    def get_case_files_by_data_category(self, project, datcat):
        basename = "%s_%s" %(project, datcat.replace(" ", "-"))
        odir = os.path.join(self.out, "%s" %(basename) )
        if( not os.path.exists(odir) ):
            os.makedirs(odir)

        fsodir = os.path.join(self.out, "%s" %(basename), "files" )
        if( not os.path.exists(fsodir) ):
            os.makedirs(fsodir)

        opath = os.path.join(odir, "files_metadata.tsv")
        if( not os.path.isfile(opath) ):
            formatt = self.formats_datCategory[datcat]
            
            filters = { "op": "and",
                    "content":[
                        {
                            "op":"in",
                            "content":{
                                "field":"cases.project.project_id",
                                    "value": [ project ]
                            }
                        },
                        {
                            "op":"=",
                            "content":{
                                "field":"data_category",
                                    "value": datcat
                            }
                        },
                        {
                            "op":"=",
                            "content":{
                                "field":"access",
                                    "value": "open"
                            }
                        },
                        {
                            "op":"in",
                            "content":{
                                "field":"data_format",
                                    "value": formatt
                            }
                        }
                    ]
                }

            fields = [
                "file_id",
                "file_size",
                "platform",
                "cases.project.project_id",
                "cases.submitter_id",
                "cases.case_id",
                "cases.samples.tumor_descriptor",
                "cases.samples.tissue_type",
                "cases.demographic.ethnicity",
                "cases.demographic.gender",
                "cases.demographic.race"
                "cases.diagnoses.age_at_diagnosis",
                "cases.project.primary_site",
                "cases.project.disease_type"
            ]

            fields = ",".join(fields)
            endpoint = "https://api.gdc.cancer.gov/files"

            params = {
                "filters": json.dumps(filters),
                "fields": fields,
                "format": "TSV",
                "size": "2000"
            }

            response = requests.get(endpoint, params = params)
            #opath = 'tmp.tsv'
            f = open( opath, 'w')
            f.write(response.text)
            f.close()

        file_list = set()
        df = pd.read_csv(opath, sep="\t")
        file_list = df.file_id.unique()

        return odir, fsodir, file_list

    def _get_file_by_uuid(self, fsodir, uuid):
        opath = os.path.join( fsodir, "raw_%s.out" %(uuid) )
        if( not os.path.isfile(opath) ):
            url = "https://api.gdc.cancer.gov/data/%s" %(uuid)
            '''
            response = requests.get(url)
            f = open( opath, 'w')
            f.write(response.text)
            f.close()
            '''
            subprocess.run(["wget", "-O", opath, url])

            try:
                with gzip.open(opath, 'rt') as f:
                    dat = f.readlines()
            except:
                with open(opath, 'r') as f:
                    dat = f.readlines()
                    
            f = open( opath,'w')
            f.write(''.join(dat))
            f.close()

    def _get_clinical_properties(self, path):
        uuid = path.split('/')[-1].split('.')[0].replace('raw_','')
        dat = { 'uuid': uuid }
        tree = ET.parse(path)
        root = tree.getroot()
        nodes = list(root.iter())
        
        targets = ['tumor_tissue_site', 'pathologic_stage', 'histological_type', 'days_to_death', 'age_at_initial_pathologic_diagnosis', 'days_to_last_known_alive', 'days_to_last_followup', 'days_to_last_follow_up', 'gender', 'ethnicity', 'race', 'vital_status']
        for t in targets:
            value = list( filter(lambda x: x.tag.endswith(t), nodes) )
            try:
                dat[t] = value[0].text
            except:
                dat[t] = ''

        dat['qty_follow_ups'] = len( list( filter(lambda x: x.tag.endswith('follow_up'), nodes) ) )

        # Information about treatment based on drugs
        drugs = list( filter(lambda x: x.tag.endswith('drugs'), nodes) )
        delements = list(drugs[0].iter())
        treat_drugs = []
        for d in delements:
            if( d.tag.endswith('drug') ):
                props = list( d.iter() )
                #print(props)
                try:
                    dname = list( filter(lambda x: x.tag.endswith('drug_name'), props) )[0].text
                except:
                    dname = ''

                try:
                    dtype = list( filter(lambda x: x.tag.endswith('therapy_type'), props) )[0].text
                except:
                    dtype = ''

                try:
                    start = list( filter(lambda x: x.tag.endswith('days_to_drug_therapy_start'), props) )[0].text
                except:
                    start = ''

                try:
                    end = list( filter(lambda x: x.tag.endswith('days_to_drug_therapy_end'), props) )[0].text
                except:
                    end = ''

                try:
                    duration = int(end) - int(start)
                except:
                    duration = ''
                
                try:
                    dosis_value = list( filter(lambda x: x.tag.endswith('prescribed_dose'), props) )[0].text
                except:
                    dosis_value = ''
                
                try:
                    dosis_unit = list( filter(lambda x: x.tag.endswith('prescribed_dose_units'), props) )[0].text
                except:
                    dosis_unit = ''
                obj = { 'type': dtype, 'name': dname, 'start': start, 'end': end, 'duration': duration, 'dosis': dosis_value, 'dosis_unit': dosis_unit }
                #print(obj)
                treat_drugs.append( obj )

        dat['qty_drugs'] = len(treat_drugs)
        dat['drug_details'] = treat_drugs

        # Information about treatment based on radiation
        rads = list( filter(lambda x: x.tag.endswith('radiations'), nodes) )
        delements = list(rads[0].iter())
        treat_rads = []
        for d in delements:
            if( d.tag.endswith('radiation') ):
                props = list( d.iter() )
                #print(props)
                try:
                    dtype = list( filter(lambda x: x.tag.endswith('radiation_type'), props) )[0].text
                except:
                    dtype = ''

                try:
                    start = list( filter(lambda x: x.tag.endswith('days_to_radiation_therapy_start'), props) )[0].text
                except:
                    start = ''

                try:
                    end = list( filter(lambda x: x.tag.endswith('days_to_radiation_therapy_end'), props) )[0].text
                except:
                    end = ''

                try:
                    duration = int(end) - int(start)
                except:
                    duration = ''
                
                try:
                    dosis_value = list( filter(lambda x: x.tag.endswith('radiation_dosage'), props) )[0].text
                except:
                    dosis_value = ''
                
                try:
                    dosis_unit = list( filter(lambda x: x.tag.endswith('units'), props) )[0].text
                except:
                    dosis_unit = ''

                obj = { 'type': dtype, 'start': start, 'end': end, 'duration': duration, 'dosis': dosis_value, 'dosis_unit': dosis_unit }
                #print(obj)
                treat_rads.append( obj )

        dat['qty_radiation'] = len(treat_rads)
        dat['radiation_details'] = treat_rads

        return dat

    def extract_data_clinical(self, odir, fsodir):
        '''
        opath = os.path.join(odir, "main_table.tsv")
        header = ["sample_id", "race", "ethnicity", "gender","tumor_tissue_site", "pathologic_stage", "days_to_death", "age_at_diagnosis", "qty_follow_ups", "tumor_status", "vital_status"]

        lines = [ header ]
        f = open(opath, 'w')
        lines = list( map( lambda x: '\t'.join(x), lines ))
        lines = '\n'.join(lines)
        f.write(lines+'\n')
        f.close()
        '''

        opath = os.path.join(odir, "data_cases.json")
        if( not os.path.exists(opath) ):
            dat_cases = []
            for f in tqdm(os.listdir(fsodir)):
                path = os.path.join(fsodir, f)
                instance = self._get_clinical_properties(path)
                dat_cases.append( instance )
            json.dump( dat_cases, open(opath, 'w') )
            #shutil.rmtree(fsodir)

    def parse_clinical_data(self, project):
        datcat = 'clinical'
        odir, fsodir, file_list = self.get_case_files_by_data_category(project, datcat)
        for uuid in file_list:
            self._get_file_by_uuid( fsodir, uuid)
        self.extract_data_clinical(odir, fsodir)

    def _get_snprs_annotation(self, path):
        dat = {}
        feat_cols = ['MAX_AF', 'MAX_AF_POPS', 'DOMAINS', 'Hugo_Symbol', 'SWISSPROT', 'Variant_Classification', 'Consequence', 'IMPACT', 'VARIANT_CLASS', 'dbSNP_RS', 'SIFT', 'PolyPhen', 'CLIN_SIG']
        df = pd.read_csv(path, sep='\t', comment='#')
        df = df[ (~ df['Consequence'].str.lower().str.contains('synonymous')) ]
        for i in df.index:
            gene = df.loc[i, "Hugo_Symbol"]
            key = df.loc[i, "dbSNP_RS"]
            if(key == 'novel'):
                key = '%s_novel' %(gene)
            if(not key in dat):
                dat[key] = { }

            for c in feat_cols:
                cname = c.lower()
                v = str(df.loc[i, c])
                if(c == "DOMAINS"):
                    aux = v.split(';')
                    for el in aux:
                        if( el.lower().startswith('pfam') ):
                            v = el.split(':')[-1]
                dat[key][cname] = v
        del df

        return dat

    def _get_mutationSnv_properties(self, path):
        '''
        distribution combinations:
            chromosome
            Consequence
            Hugo_Symbol

        '''
        group_cols = ["Chromosome", "Hugo_Symbol", "Consequence", "locationAA"]
        uuid = path.split('/')[-1].split('.')[0].replace('raw_','')
        dat = { 'uuid': uuid, 'counts': {} }

        df = pd.read_csv(path, sep='\t', comment='#')
        #df = df[ (df['Consequence'] != 'synonymous_variant') ]
        df = df[ (~ df['Consequence'].str.lower().str.contains('synonymous')) ]
        df['locationAA'] = df['Hugo_Symbol']+'_'+df['HGVSp']

        for c in group_cols:
            keys = df[ [c, 'callers'] ].groupby(c).count().index.values
            values = df[ [c, 'callers'] ].groupby(c).count().callers.values
            normc = c.lower()
            dat['counts'][normc] = dict( zip( keys, values ) )

        del df

        return dat

    def extract_data_mutationSnv(self, odir, fsodir, project):
        mapp = self._get_map_case_file(odir)
        meta_cases = self.get_cases_info_by_project(project)

        cols_stratification = ['race','gender', 'ethnicity']

        rsdat = {}
        rspath = os.path.join(odir, 'general_snprs_info.json')
        tpath = os.path.join(odir, 'by_all_table_cases.tsv')
        if( not os.path.exists(tpath) ):
            header = ["demo_variable", "subgroup", "feature", "feature_value", "uuid", "count"] 
            for c in cols_stratification+['all']:
                k = "by_%s" %(c)
                tpath = os.path.join(odir, '%s_table_cases.tsv' %(k) )
                f = open(tpath, 'w')
                f.write( '\t'.join(header)+'\n' )
                f.close()

            for f in tqdm(os.listdir(fsodir)):
                lines = {}
                for c in cols_stratification+['all']:
                    k = "by_%s" %(c)
                    lines[k] = []

                path = os.path.join(fsodir, f)

                instance = self._get_mutationSnv_properties(path)
                uuid = instance["uuid"]

                case_id = mapp[uuid]
                if( case_id in meta_cases ):
                    rsdat[uuid] = self._get_snprs_annotation(path)
                    
                    cnts = instance["counts"]
                    metas = cnts.keys()
                    for m in metas:
                        value_cnts = cnts[m]
                        for v in value_cnts:
                            vcount = value_cnts[v].item()

                            for c in cols_stratification:
                                k = "by_%s" %(c)
                                vs = meta_cases[case_id][c]
                                lines[k].append( [k, vs, m, v, uuid, vcount] )

                            lines['by_all'].append( ['all', '-', m, v, uuid, vcount] )
                    

                    for c in cols_stratification+['all']:
                        k = "by_%s" %(c)
                        tpath = os.path.join(odir, '%s_table_cases.tsv' %(k) )
                        ls = list( map( lambda x: '\t'.join( [ str(y) for y in x ] ), lines[k] ))
                        f = open(tpath, "a")
                        f.write("\n".join(ls) + "\n")
                        f.close()
            

            json.dump( rsdat, open( rspath, 'w') )
            #shutil.rmtree(fsodir)


    def _get_map_case_file(self, odir):
        mapp = {}
        opath = os.path.join(odir, 'mapping_file_case.json')

        if( not os.path.exists(opath) ):
            path = os.path.join(odir, 'files_metadata.tsv')
            df = pd.read_csv( path, sep='\t')
            mapp = dict( zip( df['file_id'].values, df['cases.0.case_id'].values ) )
            json.dump( mapp, open(opath, 'w') )
        else:
            mapp = json.load( open(opath, 'r') )

        return mapp
    
    def _get_map_file_condition(self, odir):
        mapp = {}
        opath = os.path.join(odir, 'mapping_file_condition.json')

        if( not os.path.exists(opath) ):
            path = os.path.join(odir, 'files_metadata.tsv')
            df = pd.read_csv( path, sep='\t')
            mapp = dict( zip( df['file_id'].values, df['cases.0.samples.0.tissue_type'].values ) )
            json.dump( mapp, open(opath, 'w') )
        else:
            mapp = json.load( open(opath, 'r') )

        return mapp
    
    def parse_mutationSnv_data(self, project):
        datcat = 'simple nucleotide variation'
        odir, fsodir, file_list = self.get_case_files_by_data_category(project, datcat)
        if( len(file_list) < 600 ):
            for uuid in file_list:
                self._get_file_by_uuid( fsodir, uuid)
            self.extract_data_mutationSnv(odir, fsodir, project)
        else:
            print('Skipping big files in ', project)

    def _get_gene_counts(self, path):
        df = pd.read_csv( path, sep='\t', comment='#')
        df = df[ (~df.gene_name.isna()) ]
        cnts = dict( zip( df.gene_name.values, df.unstranded.values ))
        
        return cnts

    def extract_data_expressionCounts(self, odir, fsodir, project, accepted):
        '''
        Experimental design deseq
        counts: index-rows are the case_sample_ids and the columns are the genes
        metadata: index-rows are the case_sample_ids and the columns are informations about the samples but the unique mandatory one is the condition (normal x tumor)
        
        when using deseq, set case_sample col as index
        '''

        cols_stratification = ['race','gender', 'ethnicity']
        gene_names = set()
        
        cpath = os.path.join(odir, "deseq_table_counts.tsv")
        if( not os.path.exists(cpath) ):
            mpath = os.path.join(odir, "deseq_table_meta.tsv")
            mheader = ["case_sample", "condition"] + cols_stratification
            mlines = [ mheader ]

            mappc = self._get_map_file_condition(odir)
            meta_cases = self.get_cases_info_by_project(project)
            cnts = {}
            for case_id in accepted:
                if(case_id in meta_cases):
                    groups = []
                    for cs in cols_stratification:
                        k = "by_%s" %(cs)
                        subgroup = meta_cases[case_id][cs]
                        groups.append(subgroup)

                    for uuid in accepted[case_id]:
                        condition = mappc[uuid]
                        path = os.path.join( fsodir, "raw_%s.out" %(uuid) )
                        case_sample = "%s_%s" %(case_id, uuid)
                        mlines.append( [case_sample, condition] + groups )

                        cnts[case_sample] = self._get_gene_counts(path)
                        gene_names.update( list(cnts[case_sample]) )
            
            # writing counts
            cheader = ["case_sample"] + list(gene_names)
            clines = [ cheader ]
            for cs in cnts:
                values = [cs]
                for g in gene_names:
                    try:
                        values.append( cnts[cs][g] )
                    except:
                        values.append(0)
                clines.append(values)

            clines = list( map( lambda x: '\t'.join( [ str(y) for y in x ] ), clines ))
            f = open( cpath, "w")
            f.write("\n".join(clines) + "\n")
            f.close()

            # writing metadata
            mlines = list( map( lambda x: '\t'.join( [ str(y) for y in x ] ), mlines ))
            f = open( mpath, "w")
            f.write("\n".join(mlines) + "\n")
            f.close()

    def parse_expressionCounts_data(self):
        projects = self.select_projects_open_expressionCounts()
        
        datcat = 'transcriptome profiling'
        for p in projects:
            print('--------------',p)
            odir, fsodir, file_list = self.get_case_files_by_data_category(p, datcat)
            mapp = self._get_map_case_file(odir)

            accepted = {}
            ok_files = set()
            for uuid in file_list:
                case = mapp[uuid]
                if( case in projects[p] ):
                    self._get_file_by_uuid( fsodir, uuid)
                    if( not case in accepted ):
                        accepted[case] = []
                    accepted[case].append(uuid)
            self.extract_data_expressionCounts(odir, fsodir, p, accepted)
            
        '''
        if( len(file_list) < 600 ):
            for uuid in file_list:
                self._get_file_by_uuid( fsodir, uuid)
            #self.extract_data_expressionCounts(odir, fsodir, project)
        else:
            print('Skipping big files in ', project)
        '''

    def select_projects_open_expressionCounts(self):
        cutoff = self.cutoff_cases_expression_rnaseq
        ok_projects = {}

        projects = [ "TCGA-ACC",  "TCGA-BLCA",  "TCGA-BRCA",  "TCGA-CESC",  "TCGA-CHOL",  "TCGA-COAD",  "TCGA-DLBC",  "TCGA-ESCA",  "TCGA-GBM",  "TCGA-HNSC",  "TCGA-KICH",  "TCGA-KIRC",  "TCGA-KIRP",  "TCGA-LAML",  "TCGA-LGG",  "TCGA-LIHC",  "TCGA-LUAD",  "TCGA-LUSC",  "TCGA-MESO",  "TCGA-OV",  "TCGA-PAAD",  "TCGA-PCPG",  "TCGA-PRAD",  "TCGA-READ",  "TCGA-SARC",  "TCGA-SKCM",  "TCGA-STAD",  "TCGA-TGCT",  "TCGA-THCA",  "TCGA-THYM",  "TCGA-UCEC",  "TCGA-UCS",  "TCGA-UVM" ]

        otab = os.path.join(self.out, 'expression_overview.tsv')
        if( not os.path.exists(otab) ):
            header = ["project", "qty_all_cases", "qty_cases_with_recurrence", "qty_cases_open_with_replicates", "qty_cases_with_normal&tumor", "coverage", "case_ids"]
            lines = [ header ]
            for project in projects:
                datcat = 'transcriptome profiling'
                odir, fsodir, file_list = self.get_case_files_by_data_category(project, datcat)
                pmeta = os.path.join(odir, "files_metadata.tsv")
                df = pd.read_csv( pmeta, sep='\t')

                all_cases = len( df['cases.0.case_id'].unique() )
                cases_recurrence = len( df[ df['cases.0.samples.0.tumor_descriptor'] == 'Recurrence' ]['cases.0.case_id'].unique() )

                g = df.groupby('cases.0.case_id').count().sort_values( by='platform', ascending=False)
                cases = g[ g['platform'] > 1 ].index
                df = df[ df['cases.0.case_id'].isin(cases) ]
                ok_cases = set()
                for c in cases:
                    normal = len( df[ (df['cases.0.samples.0.tissue_type']=='Normal') & ( df['cases.0.case_id']==c ) ] )
                    tumor = len( df[ (df['cases.0.samples.0.tissue_type']=='Tumor') & ( df['cases.0.case_id']==c ) ] )
                    if( normal > 0 and tumor > 0 ):
                        ok_cases.add(c)
                nexpression_ok = len(ok_cases)
                cov = (nexpression_ok/all_cases)*100
                lines.append( [project, all_cases, cases_recurrence, len(cases), nexpression_ok, cov, (','.join(ok_cases)) ] )

                if( nexpression_ok > cutoff ):
                    ok_projects[project] = ok_cases

            f = open(otab, 'w')
            lines = list( map( lambda x: '\t'.join( [ str(y) for y in x ] ), lines ))
            lines = '\n'.join(lines)
            f.write(lines+'\n')
            f.close()
        else:
            df = pd.read_csv(otab, sep='\t')
            df = df[ df["qty_cases_with_normal&tumor"] > cutoff ]
            values = list( map( lambda x: x.split(","), df.case_ids.values ) )
            ok_projects = dict( zip( df.project.values, values ) )

        return ok_projects

    def test_survival_km(self, project, datcat):
        basename = "%s_%s" %(project, datcat.replace(" ", "-"))
        odir = os.path.join(self.out, "%s" %(basename) )
        path = os.path.join(odir, "data_cases.json")

        plotdir = os.path.join(odir, "survival_plots")
        if( not os.path.exists(plotdir) ):
            os.makedirs(plotdir)

        d = json.load( open(path, 'r') )
        
        demovars = ['gender', 'ethnicity', 'race']
        valids_by = { k: {  } for k in demovars }
        valids_by['all'] = { 'all': { 'times': [], 'status': [] } }

        for s in d:
            vs = s['vital_status']

            if(vs is not None):
                samp_status = False
                if( vs.lower().find('alive') != -1 ):
                    samp_status = True

                dates = []
                dt_fields = ['days_to_death', 'days_to_last_follow_up', 'days_to_last_followup', 'days_to_last_known_alive']
                for f in dt_fields:
                    try:
                        dates.append( int( s[f] ) )
                    except:
                        dates.append( 0 )

                dt = max( dates )
                if( dt != 0 ):
                    #print(vs, vs.lower().find('alive'), samp_status, dt)

                    valids_by['all']['all']['times'].append(dt)
                    valids_by['all']['all']['status'].append(samp_status)

                    for dv in demovars:
                        try:
                            v = s[dv].replace(' ','_')
                        except:
                            v = 'Not_informed'

                        if( not v in valids_by[dv] ):
                            valids_by[dv][v] = { 'times': [], 'status': [] }
                        valids_by[dv][v]['times'].append(dt)
                        valids_by[dv][v]['status'].append(samp_status)
        
        otab = os.path.join(odir, 'survival_probs.tsv')
        header = ["demographic_dimension", "subgroup", "nsamples", "mean_survival_prob", "x", "y", "ci_low", "ci_upper"]
        lines = [ header ]
        
        for dim in valids_by:
            for subgroup in valids_by[dim]:
                status = valids_by[dim][subgroup]['status']
                times = valids_by[dim][subgroup]['times']
                #print(dim, subgroup, status)

                if( len(times) > 2 ):
                    time, survival_prob, conf_int = kaplan_meier_estimator( status, times, conf_type="log-log" )
                    nsamples = len(times)
                    ms = sum(survival_prob)/len(survival_prob)
                    x = ','.join( [ str(el) for el in time ] )
                    y = ','.join( [ str(el) for el in survival_prob ] )
                    cil = ','.join( [ str(el) for el in conf_int[0] ] )
                    ciu = ','.join( [ str(el) for el in conf_int[1] ] )

                    '''
                    opath = os.path.join( plotdir, '%s_%s_km_plot.png' %(dim, subgroup) )
                    plt.step(time, survival_prob, where="post")
                    plt.fill_between(time, conf_int[0], conf_int[1], alpha=0.25, step="post")
                    plt.ylim(0, 1)
                    plt.ylabel(r"est. probability of survival $\hat{S}(t)$")
                    plt.xlabel("time $t$, nsamples = %i, mean Surv. Prob. %.2f" %(nsamples, ms) )
                    plt.savefig(opath)
                    plt.clf()
                    '''

                    ll = [ dim, subgroup, nsamples, ms, x, y, cil, ciu ]
                    lines.append( [ str(x) for x in ll ] )

                    valids_by[dim][subgroup]['nsamples'] = nsamples
                    valids_by[dim][subgroup]['mean_prob'] = ms
                    valids_by[dim][subgroup]['x'] = time.tolist()
                    valids_by[dim][subgroup]['y'] = survival_prob.tolist()
                    valids_by[dim][subgroup]['cil'] = conf_int[0].tolist()
                    valids_by[dim][subgroup]['ciu'] = conf_int[1].tolist()

        f = open(otab, 'w')
        lines = list( map( lambda x: '\t'.join(x), lines ) )
        lines = '\n'.join(lines)
        f.write(lines+'\n')
        f.close()

        ojson = os.path.join(odir, 'survival_probs.json')
        json.dump(valids_by, open(ojson, 'w') )

    def _compress_files(self, project, datcat, namefile, remove=False):
        fname = namefile.split('.')[0]
        basename = "%s_%s" %(project, datcat.replace(" ", "-"))

        dest = '/aloy/home/ymartins/data_processed/' # change to self.out
        dest = self.out
        odir = os.path.join( dest, "%s" %(basename) )

        cwd = os.getcwd()
        fsodir = os.path.join( dest, "%s" %(basename), namefile )
        if( os.path.exists(fsodir) ):
            os.chdir(odir)
            os.system("tar -zcf %s.tar.gz %s" %(fname, namefile) )
            os.chdir(cwd)

            if(remove):
                shutil.rmtree(fsodir)

    def _decompress_files(self, project, datcat, namefile):
        fname = namefile.split('.')[0]
        basename = "%s_%s" %(project, datcat.replace(" ", "-"))
        odir = os.path.join(self.out, "%s" %(basename) )

        cwd = os.getwd()
        fsodir = os.path.join(self.out, "%s" %(basename), namefile )
        os.chdir(odir)
        os.system("tar -zxf %s.tar.gz" %(fname) )
        os.chdir(cwd)

    def run(self):
        '''
        option = sys.argv[1]
        print(option)
        #try:
        eval('self.run_'+option)()
        #except:
        #    print("Function not available")
        '''
        p = 'TCGA-BRCA'
        #self.get_vital_status_files(p)

        p = 'TCGA-ACC'
        dc = 'clinical'
        #self.get_cases_info_by_project(p)
        #self.parse_mutationSnv_data(p)
        #self.parse_expressionCounts_data(p)

        p = 'TCGA-LIHC'
        #_ = self.get_case_files_by_data_category(p, 'transcriptome profiling')
        
        projects = [ "TCGA-ACC",  "TCGA-BLCA",  "TCGA-BRCA",  "TCGA-CESC",  "TCGA-CHOL",  "TCGA-COAD",  "TCGA-DLBC",  "TCGA-ESCA",  "TCGA-GBM",  "TCGA-HNSC",  "TCGA-KICH",  "TCGA-KIRC",  "TCGA-KIRP",  "TCGA-LAML",  "TCGA-LGG",  "TCGA-LIHC",  "TCGA-LUAD",  "TCGA-LUSC",  "TCGA-MESO",  "TCGA-OV",  "TCGA-PAAD",  "TCGA-PCPG",  "TCGA-PRAD",  "TCGA-READ",  "TCGA-SARC",  "TCGA-SKCM",  "TCGA-STAD",  "TCGA-TGCT",  "TCGA-THCA",  "TCGA-THYM",  "TCGA-UCEC",  "TCGA-UCS",  "TCGA-UVM" ]
        dc = 'transcriptome profiling'
        dc = 'simple nucleotide variation'
        for p in tqdm(projects):
            #self.parse_clinical_data(p)
            #self.test_survival_km(p, dc)
            #self._compress_files(p, dc, remove=True)
            #self._compress_files(p, dc, "table_cases.tsv", remove=False)
            
            self.parse_mutationSnv_data(p)
            
        #self.select_projects_open_expressionCounts()
        #self.parse_expressionCounts_data()

if( __name__ == "__main__" ):
    o = DataWrangler()
    o.run()
