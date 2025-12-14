import os
import sys
import json
import gzip
import shutil
import requests
import pandas as pd
from tqdm import tqdm
import xml.etree.ElementTree as ET
from matplotlib import pyplot as plt
from sksurv.nonparametric import kaplan_meier_estimator

class DataWrangler:
    def __init__(self, fout = '../data_processed'):
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
            response = requests.get(url)
            f = open( opath, 'w')
            f.write(response.text)
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

    def get_vital_status_files(self, project):

        fields = [
            "submitter_id",
            "annotations.case_id",
            "cases.annotations.status",
            "cases.diagnoses.age_at_diagnosis",
            "cases.diagnoses.days_to_death",
            "cases.diagnoses.days_to_last_follow_up",
            "cases.demographic.gender",
            "files.cases.diagnoses.vital_status",
            "cases.project.primary_site",
            "cases.project.disease_type"
            ]


        fields = ",".join(fields)


        cases_endpt = "https://api.gdc.cancer.gov/v0/cases"
        cases_endpt = "https://api.gdc.cancer.gov/files"


        filters = {
            "op": "in",
            "content":{
                "field": "cases.project.project_id",
                "value": [project]
                }
            }


        # With a GET request, the filters parameter needs to be converted
        # from a dictionary to JSON-formatted string


        params = {
            "filters": json.dumps(filters),
            "fields": fields,
            "format": "TSV",
            "size": "1000"
            }

        response = requests.get(cases_endpt, params = params)
        f=open(project+'_clinical.tsv','w')
        f.write(response.text)
        f.close()


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
        
        projects = [ "TCGA-ACC",  "TCGA-BLCA",  "TCGA-BRCA",  "TCGA-CESC",  "TCGA-CHOL",  "TCGA-COAD",  "TCGA-DLBC",  "TCGA-ESCA",  "TCGA-GBM",  "TCGA-HNSC",  "TCGA-KICH",  "TCGA-KIRC",  "TCGA-KIRP",  "TCGA-LAML",  "TCGA-LGG",  "TCGA-LIHC",  "TCGA-LUAD",  "TCGA-LUSC",  "TCGA-MESO",  "TCGA-OV",  "TCGA-PAAD",  "TCGA-PCPG",  "TCGA-PRAD",  "TCGA-READ",  "TCGA-SARC",  "TCGA-SKCM",  "TCGA-STAD",  "TCGA-TGCT",  "TCGA-THCA",  "TCGA-THYM",  "TCGA-UCEC",  "TCGA-UCS",  "TCGA-UVM" ]
        for p in tqdm(projects):
        
            self.parse_clinical_data(p)
            self.test_survival_km(p, dc)

if( __name__ == "__main__" ):
    o = DataWrangler()
    o.run()
