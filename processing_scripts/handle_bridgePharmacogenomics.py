import os
import json
import gzip
import requests

import numpy as np
import pandas as pd
import gseapy as gp

from tqdm import tqdm
from time import sleep
from scipy import stats

from process_data import DataWrangler

"""
Steps inside a project:
get cases drugs
get cases mutations
    raw cols: Hugo_Symbol, HGVSp_Short, Consequence
get civic molecular profiles
    raw cols: molecular_profile (split by ' OR '), therapies (split by comma), evidence_type, evidence_direction, significance
    put relation mutation -> response/resistance -> drug, and the cancer site
"""

class HandleBridgePharcogenomicsAnalysis:
    def __init__(self, fout):
        self.data_dir = "/var/www/html/gdc_explorer_app/data_processed/"
        fout = self.data_dir

        self.projects = list( filter( lambda x: x.startswith('TCGA'), os.listdir( self.data_dir ) ))

        self.out = fout
        if( not os.path.exists(self.out) ):
            os.makedirs( fout )

        self.proc = DataWrangler(fout)
        self.map_disease, self.revmap_disease = self._get_map_project_diseaseCivic()

    def _get_map_project_diseaseCivic(self):
        projects = [ "TCGA-ACC",  "TCGA-BLCA",  "TCGA-BRCA",  "TCGA-CESC",  "TCGA-CHOL",  "TCGA-COAD",  "TCGA-DLBC",  "TCGA-ESCA",  "TCGA-GBM",  "TCGA-HNSC",  "TCGA-KICH",  "TCGA-KIRC",  "TCGA-KIRP",  "TCGA-LAML",  "TCGA-LGG",  "TCGA-LIHC",  "TCGA-LUAD",  "TCGA-LUSC",  "TCGA-MESO",  "TCGA-OV",  "TCGA-PAAD",  "TCGA-PCPG",  "TCGA-PRAD",  "TCGA-READ",  "TCGA-SARC",  "TCGA-SKCM",  "TCGA-STAD",  "TCGA-TGCT",  "TCGA-THCA",  "TCGA-THYM",  "TCGA-UCEC",  "TCGA-UCS",  "TCGA-UVM" ]
        diseases = [ ["Adrenocortical Carcinoma"], ["Bladder Carcinoma"], ["Breast Cancer"], ["Cervical Cancer"], ["Cholangiocarcinoma"], ["Colon Adenocarcinoma", "Colorectal Cancer"], ["Diffuse Large B-cell Lymphoma"], ["Esophageal Carcinoma", "Esophagus Squamous Cell Carcinoma"], ["Glioblastoma"], ["Head And Neck Squamous Cell Carcinoma"], ["Chromophobe Renal Cell Carcinoma", "Renal cell carcinoma"], ["Kidney Clear Cell Sarcoma", "Renal cell carcinoma"], ["Papillary Renal Cell Carcinoma", "Renal cell carcinoma"], ["Acute Myeloid Leukemia"], ["Low Grade Glioma"], ["Hepatocellular Carcinoma"], ["Lung Adenocarcinoma"], ["Lung Squamous Cell Carcinoma"], ["Malignant Pleural Mesothelioma"], ["Ovarian Cancer"], ["Pancreatic Cancer"], ["Pheochromocytoma", "Paraganglioma"], ["Prostate Carcinoma"], ["Rectum Cancer"], ["Sarcoma"], ["Skin Melanoma"], ["Stomach Cancer"], ["Testicular Cancer"], ["Thyroid Cancer"], ["Thymic Carcinoma"], ["Uterine Corpus Endometrial Carcinoma", "Endometrial Cancer"], ["Uterine Cancer"], ["Uveal Melanoma"] ]
        mapp = dict( zip( projects, diseases ))
        
        revmapp = {}
        for k in mapp:
            arr = mapp[k]
            for d in arr:
                if( not d in revmapp):
                    revmapp[d] = []
                revmapp[d].append(k)

        return mapp, revmapp

    def get_civic_data(self):
        dat = {}

        dir = '../external_db'
        path = os.path.join(dir, 'cividb_jun-26.tsv')
        df = pd.read_csv(path, sep='\t', comment='#')
        df = df[ (df["evidence_type"] == 'Predictive') & (df['evidence_direction']=='Supports') ]

        for i in df.index:
            mp = df.loc[i, 'molecular_profile']
            disease = df.loc[i, 'disease']
            mps = list(map( lambda x: x.replace(' ', '-'), mp.split(' OR ')))
            drugs = df.loc[i, 'therapies']
            if( str(drugs) != 'nan' ):
                drugs = drugs.split(',')
                signal = df.loc[i, 'significance']

                for _id in mps:
                    if( not _id in dat ):
                        dat[_id] = []
                    for d in drugs:
                        dat[_id].append( { 'drug': d, 'significance': signal, 'disease': disease } )
            else:
                print(df.loc[i, :].values)

        return dat

    def _get_snprs_annotation(self, path):
        dat = {}
        feat_cols = ['Reference_Allele', 'Tumor_Seq_Allele1', 'Tumor_Seq_Allele2', 'HGVSc', 'MAX_AF', 'MAX_AF_POPS', 'DOMAINS', 'Hugo_Symbol', 'SWISSPROT', 'Variant_Classification', 'Consequence', 'IMPACT', 'VARIANT_CLASS', 'dbSNP_RS', 'SIFT', 'PolyPhen', 'CLIN_SIG', "Chromosome"]
        df = pd.read_csv(path, sep='\t', comment='#')
        df = df[ (~ df['Consequence'].str.lower().str.contains('synonymous')) ]
        df['locationAA'] = df['Hugo_Symbol']+'-'+df['HGVSp_Short']
        df = df[ ~df['locationAA'].isna() ]
        df['locationAA'] = df['locationAA'].map( lambda x: x.replace('p.','') )
        for i in df.index:
            aachange = df.loc[i, 'HGVSp']
            gene = df.loc[i, "Hugo_Symbol"]
            rsid = df.loc[i, "dbSNP_RS"]
            t_ref = df.loc[i, "t_ref_count"]
            t_alt = df.loc[i, "t_alt_count"]
            if(rsid not in ['novel', 'nan'] ):
                #rsid = '%s_novel' %(gene)

                key = df.loc[i, "locationAA"]
                if( str(key) != 'nan' ):
                    if(not key in dat):
                        dat[key] = { }

                    dat[key]["aa_change"] = aachange
                    for c in feat_cols:
                        cname = c.lower()
                        v = str(df.loc[i, c]).split('(')[0]

                        if(c == "DOMAINS"):
                            aux = v.split(';')
                            for el in aux:
                                if( el.lower().startswith('pfam') ):
                                    v = el.split(':')[-1]
                        dat[key][cname] = v
        del df

        return dat

    def merge_datasources(self):
        feat_cols = ['Variant_Classification', 'Consequence', 'PolyPhen', 'CLIN_SIG', "Chromosome"]
        feat_cols = list(map( lambda x: x.lower(), feat_cols ))
        header = ["case_id", "tissue_type", "civic_gene_hits", "civic_mutation_hits", "civic_drug_hits", "civic_mut&drug_hits", "civic_matched_drugs_for_mutation", "civic_matched_diseases_for_mutation", "gene", "mutation", "drug", "cvdb_significance", "cvdb_disease", "race", "gender", "ethnicity"] + feat_cols

        cvdb = self.get_civic_data()
        cvdb_muts = list(cvdb)
        cvdb_genes = set( map( lambda x: x.split('-')[0], cvdb_muts ))
        cvdb_drugs = set()
        for m in cvdb_muts:
            cvdb_drugs.update( list(map( lambda x: x['drug'], cvdb[m] )) )

        for project in tqdm(self.projects):
            indir = os.path.join(self.data_dir, project)
            dcases = self.proc._get_cases_metadata(project)

            datcat = 'simple nucleotide variation'
            odir, fsodir, file_list = self.proc.get_case_files_by_data_category(project, datcat)
            lines = [header]

            mapp = self.proc._get_map_case_file(odir)
            mapp_tissue = self.proc._get_map_file_condition(odir)
            for uuid in file_list:
                f = "raw_%s.out" %(uuid)
                path = os.path.join(fsodir, f)
                case_id = mapp[uuid]
                tissue_type = mapp_tissue[uuid]

                if(case_id in dcases):
                    dclin = dcases[case_id]
                    snps = self._get_snprs_annotation(path)

                    drg_details = dclin["drug_details"]
                    drugs = set()
                    for d in drg_details:
                        if('name' in d):
                            drugs.add(d['name'])

                    gender = dclin["gender"]
                    race = dclin["race"]
                    ethnicity = dclin["ethnicity"]
                    muts = list(snps)
                    for m in muts:
                        gene = m.split('-')[0]
                        flag_cvdb_genes = (gene in cvdb_genes)
                        flag_cvdb_muts = (m in cvdb_muts)

                        datmuts = []
                        for f in feat_cols:
                            datmuts.append( snps[m][f] )

                        for dr in drugs:
                            flag_cvdb_drugs = (dr in cvdb_drugs)
                            flag_cvdb_both = (flag_cvdb_muts and flag_cvdb_drugs)
                            
                            signal = '-'
                            disease = '-'
                            matched_drugs = '-'
                            matched_dis = '-'
                            if( m in cvdb):
                                matched_drugs = list( map( lambda x: x['drug'], cvdb[m] ))
                                matched_drugs = ','.join(matched_drugs)
                                matched_dis = list( map( lambda x: x['disease'], cvdb[m] ))
                                matched_dis = ','.join(matched_dis)
                                fcv = list(filter( lambda x: x['drug']==dr, cvdb[m] ))
                                if( len(fcv) > 0 ):
                                    signal = fcv[0]['significance']
                                    disease = fcv[0]['disease']

                            el = [case_id, tissue_type, flag_cvdb_genes, flag_cvdb_muts, flag_cvdb_drugs, flag_cvdb_both, matched_drugs, matched_dis, gene, m, dr, signal, disease, race, gender, ethnicity] + datmuts
                            lines.append(el)
                else:
                    print('not found case --> ', project, case_id)

            opath = os.path.join(odir, "data_merge_pharmaco_mutations.tsv")
            lines = list( map( lambda x: '\t'.join( [ str(y) for y in x ] ), lines ))
            f = open( opath, "w")
            f.write("\n".join(lines) + "\n")
            f.close()

    def generate_app_pgx_compiled_data(self):
        '''
        Interface planning:
        - a plot with the mutations in the cases matched with civicdb, x=mutation, y=cases count
            A button below showing in a table the other mutations not shown in plot for reasons of clarity

        when clicking on a mutation:
            - a table with relevant information for each case containing the mutation, prescribed drugs x drugs associated with mutation in civic and what the mutation cause (response or resistance) |
            - each line of the table there is a button to open specific analysis area for the cases, and shows the overview of mutation features found for the case
            
            ok- plots and info for general overview: distribution of cases by gender, race, all civic diseases | drug-significance pairs.
            
            - columns of cases table: case_id, demo info ( race| gender | ethnicity), chromosomes, genes, prescribed drugs in treatment, number of muts in civic & all muts
        
        when clicking in details button in a specific case:
            - plot overview of cases-specific information: distribution of genes, chromosome, consequence, drugs used in treatment
            - button for the user to open the full mutation list in table for the case, and flag the ones that are in civic

        maybe as an accordeon, could have a section showing the drugs plot, that are in the cases and civic db
        when clicking in a drug, it also appear the same cases table.
        If mapping from name to atc code, it is possible to get data of pharmoGenomics from pgxdb ( https://pgx-db.org/rest-api/atc/pgx/D07XB05/ ) - important properties: Phenotype_Category (efficacy), Variant_or_Haplotypes (rsid of mutation), Direction_of_effect (increased), PD_PK_terms (response to), PMID

        --- interface
        table pagination = https://github.com/wenzhixin/bootstrap-pagination
        modal = https://getbootstrap.com/docs/5.3/components/modal/

        '''
        cvdb = self.get_civic_data()

        for project in tqdm(self.projects):
            indir = os.path.join(self.data_dir, project)
            dcases = self.proc._get_cases_metadata(project)

            datcat = 'simple nucleotide variation'
            odir, fsodir, file_list = self.proc.get_case_files_by_data_category(project, datcat)
            path = os.path.join(odir, "data_merge_pharmaco_mutations.tsv")
            df = pd.read_csv( path, sep='\t')
            
            #df[ df['civic_mutation_hits'] ].case_id.unique() # number of cases that had a mutation annotated in civic
            #df[ df['civic_mutation_hits'] ].groupby(['race', 'gender', 'ethnicity']).count() # same filter but grouping the table by race, gender and ethnicity

            # --- make mapp of case to all mutations observed, in or not in civic
            aux = {}
            for c in df.case_id.unique():
                aux[c] = df[ df.case_id == c ].mutation.unique()

            df = df[ df['civic_mutation_hits'] ]
            mutdata = dict() # Textual tag p: chromosome, gene, consequence, polyphen, clinsig, civic diseases, drug-variant_association
            mutcount = dict()
            casedata = dict()
            tissue_site = ""
            hist_type = ""
            for i in df.index:
                mut = df.loc[i, 'mutation']
                case = df.loc[i, 'case_id']
                tissue_type = df.loc[i, 'tissue_type']

                # filling mut info json
                if(not mut in mutdata):
                    mutcount[mut] = { "cases": 0, "race": {}, "gender": {}, "stage": {}, "tissue_type": {} }
                    mutdata[mut] = { "race": {}, "gender": {}, "stage": {}, "tissue_type": {}, "civic_diseases": {}, "drug-variant_association": {}, "cases": set() }

                mutdata[mut]["cases"].add(case)
                mutcount[mut]["cases"] = len( mutdata[mut]["cases"] )

                race = df.loc[i, 'race']
                if(not race in mutdata[mut]["race"] ):
                    mutdata[mut]["race"][race] = set()
                mutdata[mut]["race"][race].add(case)
                mutcount[mut]["race"][race] = len( mutdata[mut]["race"][race] )

                gender = df.loc[i, 'gender']
                if(not gender in mutdata[mut]["gender"] ):
                    mutdata[mut]["gender"][gender] = set()
                mutdata[mut]["gender"][gender].add(case)
                mutcount[mut]["gender"][gender] = len( mutdata[mut]["gender"][gender] )

                chromosome = str(df.loc[i, 'chromosome'])
                consequence = str(df.loc[i, 'consequence'])
                polyphen = str(df.loc[i, 'polyphen'])
                clinsig = str(df.loc[i, 'clin_sig'])
                mutdata[mut]["chromosome"] = chromosome
                mutdata[mut]["consequence"] = consequence
                mutdata[mut]["polyphen"] = polyphen
                mutdata[mut]["clinsig"] = clinsig

                cvdat = cvdb[mut]
                diseases = set( map( lambda x: x["disease"], cvdat ))
                combs = set( map( lambda x: x["significance"]+" to "+x["drug"], cvdat ))
                mutdata[mut]["civic_diseases"] = ', '.join(diseases)
                mutdata[mut]["drug-variant_association"] = ', '.join(combs)

                # filling case info json
                gene = mut.split('-')[0]
                mutdata[mut]["gene"] = gene
                
                drug = str(df.loc[i, 'drug'])
                dclin = dcases[case]
                stage = dclin["pathologic_stage"]
                tissue_site = dclin["tumor_tissue_site"]
                hist_type = dclin["histological_type"]
                if(not case in casedata):
                    drg_details = dclin["drug_details"]
                    drugs = set()
                    for d in drg_details:
                        if('name' in d):
                            drugs.add(d['name'])
                    gender = dclin["gender"]
                    race = dclin["race"]
                    ethnicity = dclin["ethnicity"]

                    casedata[case] = { "race": race, "gender": gender, "ethnicity": ethnicity, "stage": stage, "prescribed_drugs": set(), "genes": set(), "muts": set(), "chromosomes": set() }
                casedata[case]["chromosomes"].add(chromosome)

                if(not stage in mutdata[mut]["stage"] ):
                    mutdata[mut]["stage"][stage] = set()
                mutdata[mut]["stage"][stage].add(case)
                mutcount[mut]["stage"][stage] = len( mutdata[mut]["stage"][stage] )

                if(not tissue_type in mutdata[mut]["tissue_type"] ):
                    mutdata[mut]["tissue_type"][tissue_type] = set()
                mutdata[mut]["tissue_type"][tissue_type].add(case)
                mutcount[mut]["tissue_type"][tissue_type] = len( mutdata[mut]["tissue_type"][tissue_type] )

                flagg = df.loc[i, 'civic_gene_hits']
                if(flagg):
                    gene += "*"
                casedata[case]["genes"].add(gene)
                
                flagd = df.loc[i, 'civic_drug_hits']
                if(flagd):
                    drug += "*"
                casedata[case]["prescribed_drugs"].add(drug)
                
                casedata[case]["muts"].add(mut)
                
            tmp = {}
            for k in mutdata:
                tmp[k] = {}
                tmp[k]["cases"] = list(mutdata[k]["cases"])
                for c in mutdata[k]:
                    if(c not in ["race", "gender", "cases", "stage", "tissue_type"] ):
                        tmp[k][c] = mutdata[k][c]
            mutdata = tmp
            
            for k in casedata:
                casedata[k]["chromosomes"] = list(casedata[k]["chromosomes"])
                casedata[k]["genes"] = list(casedata[k]["genes"])
                casedata[k]["muts"] = list(casedata[k]["muts"])
                casedata[k]["allmuts"] = list(aux[k])
                casedata[k]["prescribed_drugs"] = list(casedata[k]["prescribed_drugs"])

            print(project, mutdata, mutcount, casedata)
            alldat = { "tissue_site": tissue_site, "hist_type": hist_type, "mutation_details": mutdata, "mutation_count": mutcount, "cases_data": casedata }
            opath = os.path.join(odir, "pgx_data.json")
            json.dump( alldat, open(opath, 'w') )
            if( os.path.getsize(opath) == 2 ):
                print(project)

    def _prepare_drugs_list(self):
        drugs = set()
        cvdb = self.get_civic_data()
        for mut in cvdb:
            cvdat = cvdb[mut]
            drugs.update( list( map( lambda x: x["drug"], cvdat )) )

        for project in tqdm(self.projects):
            indir = os.path.join(self.data_dir, project)
            dcases = self.proc._get_cases_metadata(project)

            datcat = 'simple nucleotide variation'
            odir, fsodir, file_list = self.proc.get_case_files_by_data_category(project, datcat)
            path = os.path.join(odir, "data_merge_pharmaco_mutations.tsv")
            df = pd.read_csv( path, sep='\t')
            drugs.update( df.drug.values )

        drugs = list( filter( lambda x: type(x) == str, drugs ) )
        opath = os.path.join( self.out, "all_drugs.txt")
        f = open( opath, 'w')
        f.write( '\n'.join(drugs) )
        f.close()

    def get_smiles_from_graphql_pdb(self):
        """
        the file /mnt/yasdata/home/yasmmin/Dropbox/portfolio_2025/gdc_explorer_app_web/gdc_exploration_project/map_ligands_smiles.tsv has the mapping between ligand id pdb to smiles.
        To get all available durgbank ligand entries in pdb, go to advanced search, chemical search tab >> filters: Lineage Identifier - ATC (WHO) is L ; Resource Name is Drugbank
        The antibody-based drugs such as bevacizumab cannot be found this way, only the inhibitors, because they are not ligands or small molecules

        Ligands 3-letter ids:
        032,07J,09L,0LI,0WM,0XZ,1C9,1E8,1K5,1LT,1N1,1XJ,2HB,2K2,2TA,2YQ,3E8,3EW,3JD,3JW,3WF,40L,4BM,4J8,4MK,5AE,5OG,5P8,5SF,62G,69Q,6E2,6GY,6K9,6T3,6V8,6ZV,6ZZ,72Q,78P,7GI,8JC,8ZF,9CR,9JI,9NQ,9RA,9TP,9UO,A1AAC,A1B7W,A1CI6,A1D8L,A1EK0,A1ESZ,A1JCR,A1JRF,A1JS8,A1JSO,A26,A4I,A9L,AER,AQ4,AR3,ASW,AV9,AXI,AY7,B49,BAX,BBJ,BLM,BO2,C6F,C87,CBL,CEL,CFB,CL9,CP0,CPT,CTX,D16,DB8,DCF,DES,DM1,DM2,DM5,DS9,DX4,ECT,EFD,EMH,EOU,EUI,EVP,EXM,FK5,FMM,FVT,G65,GEO,GZX,HFT,HMT,HSM,I0V,ICQ,IRE,IV3,J33,J8C,JEU,JGQ,LBH,LBM,LEV,LON,LQQ,LYA,MGB,MI1,MIX,MOA,MTX,NHY,NIL,O6U,P06,P30,P31,P7D,PM6,PWV,Q4J,QO7,QOM,QWP,R1Q,RAP,REA,RPB,RQ3,RXT,S4R,SHH,STI,T0R,TA1,TIZ,TTC,TXL,TZ0,UKI,URF,VGH,VH6,VIS,VLB,VUC,XIN,XQQ,Y7W,YMX,YY3,ZD6,ZIY
        https://data.rcsb.org/graphiql/index.html?query=query%20molecule%20(%24id%3A%20String!)%20%7B%0A%20%20%20%20chem_comp(comp_id%3A%24id)%7B%0A%20%20%20%20%20%20%20%20chem_comp%20%7B%0A%20%20%20%20%20%20%20%20%20%20%20%20id%0A%20%20%20%20%20%20%20%20%20%20%20%20name%0A%20%20%20%20%20%20%20%20%20%20%20%20formula%0A%20%20%20%20%20%20%20%20%20%20%20%20pdbx_formal_charge%0A%20%20%20%20%20%20%20%20%20%20%20%20formula_weight%0A%20%20%20%20%20%20%20%20%20%20%20%20type%0A%20%20%20%20%20%20%20%20%7D%0A%20%20%20%20%20%20%20%20pdbx_reference_molecule%20%7B%0A%20%20%20%20%20%20%20%20%20%20%20%20prd_id%0A%20%20%20%20%20%20%20%20%20%20%20%20chem_comp_id%0A%20%20%20%20%20%20%20%20%20%20%20%20type%0A%20%20%20%20%20%20%20%20%20%20%20%20class%0A%20%20%20%20%20%20%20%20%20%20%20%20name%0A%20%20%20%20%20%20%20%20%20%20%20%20represent_as%0A%20%20%20%20%20%20%20%20%20%20%20%20representative_PDB_id_code%0A%20%20%20%20%20%20%20%20%7D%0A%20%20%20%20%20%20%20%20rcsb_chem_comp_info%20%7B%0A%20%20%20%20%20%20%20%20%20%20%20%20atom_count%0A%20%20%20%20%20%20%20%20%20%20%20%20bond_count%0A%20%20%20%20%20%20%20%20%20%20%20%20bond_count_aromatic%0A%20%20%20%20%20%20%20%20%20%20%20%20atom_count_chiral%0A%20%20%20%20%20%20%20%20%20%20%20%20initial_deposition_date%0A%20%20%20%20%20%20%20%20%20%20%20%20revision_date%0A%20%20%20%20%20%20%20%20%7D%0A%20%20%20%20%20%20%20%20rcsb_chem_comp_descriptor%20%7B%0A%20%20%20%20%20%20%20%20%20%20%20%20InChI%0A%20%20%20%20%20%20%20%20%20%20%20%20InChIKey%0A%20%20%20%20%20%20%20%20%20%20%20%20SMILES%0A%20%20%20%20%20%20%20%20%20%20%20%20SMILES_stereo%0A%20%20%20%20%20%20%20%20%7D%0A%20%20%20%20%20%20%20%20pdbx_reference_entity_poly_seq%20%7B%0A%20%20%20%20%20%20%20%20%20%20%20%20observed%0A%20%20%20%20%20%20%20%20%20%20%20%20mon_id%0A%20%20%20%20%20%20%20%20%20%20%20%20num%0A%20%20%20%20%20%20%20%20%7D%0A%20%20%20%20%20%20%20%20pdbx_chem_comp_identifier%20%7B%0A%20%20%20%20%20%20%20%20%20%20%20%20identifier%0A%20%20%20%20%20%20%20%20%20%20%20%20program%0A%20%20%20%20%20%20%20%20%7D%0A%20%20%20%20%20%20%20%20pdbx_chem_comp_descriptor%20%7B%0A%20%20%20%20%20%20%20%20%20%20%20%20type%0A%20%20%20%20%20%20%20%20%20%20%20%20descriptor%0A%20%20%20%20%20%20%20%20%20%20%20%20program%0A%20%20%20%20%20%20%20%20%20%20%20%20program_version%0A%20%20%20%20%20%20%20%20%7D%0A%20%20%20%20%20%20%20%20rcsb_chem_comp_synonyms%20%7B%0A%20%20%20%20%20%20%20%20%20%20name%0A%20%20%20%20%20%20%20%20%20%20type%0A%20%20%20%20%20%20%20%20%20%20provenance_source%0A%20%20%20%20%20%20%20%20%7D%20%0A%20%20%20%20%20%20%20%20drugbank%20%7B%0A%20%20%20%20%20%20%20%20%20%20%20%20drugbank_info%20%7B%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20drugbank_id%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20cas_number%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20drug_categories%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20mechanism_of_action%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20synonyms%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20name%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20drug_groups%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20description%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20affected_organisms%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20brand_names%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20indication%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20pharmacology%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20atc_codes%0A%20%20%20%20%20%20%20%20%20%20%20%20%7D%0A%20%20%20%20%20%20%20%20%20%20%20%20drugbank_target%20%7B%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20target_actions%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20name%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20interaction_type%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20seq_one_letter_code%0A%20%20%20%20%20%20%20%20%20%20%20%20%7D%0A%20%20%20%20%20%20%20%20%7D%0A%20%20%20%20%20%20%20%20rcsb_chem_comp_related%20%7B%0A%20%20%20%20%20%20%20%20%20%20%20%20resource_name%0A%20%20%20%20%20%20%20%20%20%20%20%20resource_accession_code%0A%20%20%20%20%20%20%20%20%7D%20%0A%20%20%20%20%7D%0A%7D%0A&variables=%7B%0A%20%20%22id%22%3A%20%225P8%22%0A%7D
        """
        ids = "032,07J,09L,0LI,0WM,0XZ,1C9,1E8,1K5,1LT,1N1,1XJ,2HB,2K2,2TA,2YQ,3E8,3EW,3JD,3JW,3WF,40L,4BM,4J8,4MK,5AE,5OG,5P8,5SF,62G,69Q,6E2,6GY,6K9,6T3,6V8,6ZV,6ZZ,72Q,78P,7GI,8JC,8ZF,9CR,9JI,9NQ,9RA,9TP,9UO,A1AAC,A1B7W,A1CI6,A1D8L,A1EK0,A1ESZ,A1JCR,A1JRF,A1JS8,A1JSO,A26,A4I,A9L,AER,AQ4,AR3,ASW,AV9,AXI,AY7,B49,BAX,BBJ,BLM,BO2,C6F,C87,CBL,CEL,CFB,CL9,CP0,CPT,CTX,D16,DB8,DCF,DES,DM1,DM2,DM5,DS9,DX4,ECT,EFD,EMH,EOU,EUI,EVP,EXM,FK5,FMM,FVT,G65,GEO,GZX,HFT,HMT,HSM,I0V,ICQ,IRE,IV3,J33,J8C,JEU,JGQ,LBH,LBM,LEV,LON,LQQ,LYA,MGB,MI1,MIX,MOA,MTX,NHY,NIL,O6U,P06,P30,P31,P7D,PM6,PWV,Q4J,QO7,QOM,QWP,R1Q,RAP,REA,RPB,RQ3,RXT,S4R,SHH,STI,T0R,TA1,TIZ,TTC,TXL,TZ0,UKI,URF,VGH,VH6,VIS,VLB,VUC,XIN,XQQ,Y7W,YMX,YY3,ZD6,ZIY".split(",")
        url = "https://data.rcsb.org/graphql"
        graphql_query = """

query molecule ($id: String!) {
    chem_comp(comp_id: $id){
        chem_comp {
            id
            name
            formula
            pdbx_formal_charge
            formula_weight
            type
        }
        pdbx_reference_molecule {
            prd_id
            chem_comp_id
            type
            class
            name
            represent_as
            representative_PDB_id_code
        }
        rcsb_chem_comp_info {
            atom_count
            bond_count
            bond_count_aromatic
            atom_count_chiral
            initial_deposition_date
            revision_date
        }
        rcsb_chem_comp_descriptor {
            InChI
            InChIKey
            SMILES
            SMILES_stereo
        }
        pdbx_reference_entity_poly_seq {
            observed
            mon_id
            num
        }
        pdbx_chem_comp_identifier {
            identifier
            program
        }
        pdbx_chem_comp_descriptor {
            type
            descriptor
            program
            program_version
        }
        rcsb_chem_comp_synonyms {
          name
          type
          provenance_source
        } 
        drugbank {
            drugbank_info {
                drugbank_id
                cas_number
                drug_categories
                mechanism_of_action
                synonyms
                name
                drug_groups
                description
                affected_organisms
                brand_names
                indication
                pharmacology
                atc_codes
            }
            drugbank_target {
                target_actions
                name
                interaction_type
                seq_one_letter_code
            }
        }
        rcsb_chem_comp_related {
            resource_name
            resource_accession_code
        } 
    }
}

        """
        opath = os.path.join(self.out, "map_drugbank_ligPdb_smiles.tsv")
        if( not os.path.exists(opath) ):
            header = ["lig", "name", "smiles"]
            lines = [ header ]
            for _id in tqdm(ids):
                variables = {"id": _id}
                payload = {
                    "query": graphql_query,
                    "variables": variables
                }
                response = requests.post( url, data = json.dumps(payload), headers = {'Content-Type': 'application/json'} )
                
                # Check and display the parsed JSON data
                if response.status_code == 200:
                    data = response.json()["data"]
                    lig = data["chem_comp"]["chem_comp"]["id"]
                    name = data["chem_comp"]["drugbank"]["drugbank_info"]["name"]
                    smiles = data["chem_comp"]["rcsb_chem_comp_descriptor"]["SMILES"]
                    lines.append([lig, name, smiles])
                else:
                    print(f"Query failed with status code {response.status_code}")

            lines = list( map( lambda x: '\t'.join( [ str(y) for y in x ] ), lines ))
            f = open( opath, "w")
            f.write("\n".join(lines) + "\n")
            f.close()

        df = pd.read_csv(opath, sep='\t')
        smis = df.smiles.values
        opath = os.path.join( self.out, "smiles_for_admetpred.txt")
        f = open( opath, 'w')
        f.write( '\n'.join(smis) )
        f.close()

    def get_smiles_from_chembl(self):
        path = os.path.join( self.out, "all_drugs.txt")
        drugs = open( path, 'r' ).read().split('\n')

        import sqlite3
        db_path = '/mnt/yasdata/home/yasmmin/Dropbox/portfolio_2025/gdc_explorer_app_web/gdc_exploration_project/Unconfirmed 402173.crdownload' 
        conn = sqlite3.connect(db_path)
        cursor = conn.cursor()

        # 2. Define the target compound name (e.g., 'Aspirin' or 'Ibuprofen')
        compound_name = 'Aspirin'

        # 3. Execute the SQL query
        query = """
        SELECT md.molregno, md.pref_name, cs.canonical_smiles
        FROM molecule_dictionary md
        JOIN compound_structures cs ON md.molregno = cs.molregno
        WHERE md.pref_name LIKE ?
        """

        mapp = {}
        for drug in drugs:
            # We use % around the name to act as a wildcard in case of exact casing mismatches
            cursor.execute(query, (f"%{drug}%",))
            results = cursor.fetchall()

            # 4. Display the results
            for row in results:
                molregno, pref_name, smiles = row
                print(f"ChEMBL ID: {molregno}")
                print(f"Name: {pref_name}")
                print(f"SMILES: {smiles}\n")
                mapp[drug] = smiles
        opath = os.path.join( self.out, "drug_smiles.json" )

    def clean_drug_names(self):
        path = os.path.join( self.out, "all_drugs.txt")
        drugs = open( path, 'r' ).read().split('\n')
        drugs = list( map( lambda x: x.lower(), drugs ))
        ds = set()
        for d in drugs:
            d = d.lower()
            d = d.replace(' or ', ',')
            d = d.replace(' and ', ',')
            tmp = [d]
            if( ',' in d ):
                tmp = list( map( lambda x: x.strip().lower(), d.split(',') ))
            if( '+' in d ):
                tmp = list( map( lambda x: x.strip().lower(), d.split('+') ))
            if( '/' in d ):
                tmp = list( map( lambda x: x.strip().lower(), d.split('/') ))
            tmp = list( filter( lambda x: x != 'placebo', tmp ))
            ds.update(tmp)

        opath = os.path.join( self.out, "clean_all_drugs.txt")
        f = open( opath, 'w')
        f.write( '\n'.join(ds) )
        f.close()

    def get_clinpgx_drugLabels(self):
        # All clean pgx api urls = https://api.clinpgx.org/swagger/#/

        path = os.path.join( self.out, "clean_all_drugs.txt")
        drugs = open( path, 'r' ).read().split('\n')
        dat = {}
        for d in tqdm(drugs):
            r = requests.get( f"https://api.clinpgx.org/v1/data/label?relatedChemicals.name={d}" )
            if( r.status_code == 200 ):
                dat[d] = r.json()

        opath = os.path.join(self.out, "data_all_drugs_pgx.json")
        json.dump(dat, open(opath, 'w') )

    def get_pharmvar_alleles(self):
        ids = []
        opath = os.path.join(self.out, "data_pharmvar_alleles.json")
        if( not os.path.exists(opath) ):
            headers = { "Api-Key": os.environ["pharmvar_key"] }
            r = requests.get( f"https://www.pharmvar.org/api-service/alleles/list", headers = headers )
            if( r.status_code == 200 ):
                dat = r.json()
                json.dump( dat, open(opath, 'w') )
        else:
            dat = json.load( open(opath, 'r') )
        
        return dat

    def _parse_allele_table_info(self, path):
        df = pd.read_excel(path, sheet_name="Alleles")
        
        als = list( filter( lambda x: x.startswith('*'), df.iloc[:,0].values ))
        
        i = 0
        datmeta = {}
        metas = { "rsid": "rsid", "on protein": "protein_change", "refseqgene": "refseq_genome_change", "GRCh38": "gv38_genome_change" }
        for h in df.iloc[:, 0].values:
            if( h.startswith("*") ):
                break

            hn = h.lower()
            for k in metas:
                if( hn.find(k) != -1 ):
                    hn = metas[k]

            datmeta[hn] = df.iloc[i, 1:]

            i += 1

        aldat = {}
        i = 0
        for h in df.iloc[:, 0].values:
            if( h.startswith("*") and h != "*1" ):
                vs = df.iloc[i, 1:]
                for k in datmeta:
                    j = 0
                    for v in vs:
                        if( str(v) != 'nan' ):
                            if(not h in aldat):
                                aldat[h] = {}
                            if(not k in aldat[h] ):
                                aldat[h][k] = []
                                aldat[h]["variant"] = []
                            aldat[h][k].append( datmeta[k][j] )
                            aldat[h]["variant"].append( v )
                        j+=1
            i+=1

        return aldat

    def _download_parse_allele_definition_table(self, names):
        dpath = os.path.join(self.out, "allele_definition")
        if( not os.path.isdir(dpath) ):
            os.makedirs(dpath)

        dat = {}
        opath = os.path.join(self.out, "alleles_definition.json")
        if(not os.path.exists(opath) ):
            for gn in tqdm(names):
                dat[gn] = {}
                dls = f"https://s3.pgkb.org/attachment/{gn}_allele_definition_table.xlsx"
                resp = requests.get(dls)

                opath = os.path.join( dpath, f"{gn}_alleles.xlsx")
                output = open( opath, 'wb')
                output.write(resp.content)
                output.close()

                dt = self._parse_allele_table_info(opath)
                dat[gn] = dt

            json.dump(dat, open(opath, 'w') )
        else:
            dat = json.load( open(opath, "r") )

        return dat

    def get_pharmvar_allele_info(self, ids = []):
        dpath = os.path.join(self.out, "allele_information")
        if( not os.path.isdir(dpath) ):
            os.makedirs(dpath)

        dat = {}
        opath = os.path.join( dpath, "alleles_full_information.json")
        if(not os.path.exists(opath) ):
            headers = { "Api-Key": os.environ["pharmvar_key"] }

            for i in tqdm(ids):
                popath = os.path.join(dpath, f"{i}_info.json")
                part = {}
                if( not os.path.exists(popath) ):
                    r = requests.get( f"https://www.pharmvar.org/api-service/alleles/{i}?exclude-sub-alleles=false&include-reference-variants=false&include-retired-alleles=false&include-retired-reference-sequences=false", headers = headers )
                    if( r.status_code == 200 ):
                        part = r.json()
                    else:
                        print(i, r)

                    json.dump(part, open(popath, 'w') )
                else:
                    part = json.load( open(popath, "r") )

                dat[i] = part

            json.dump(dat, open(opath, 'w') )
        else:
            dat = json.load( open(opath, "r") )

        return dat

    def process_allele_info(self):
        dat = {}

        opath = os.path.join( self.out, "allele_info_summary.json")
        if( not os.path.exists(opath) ):
            dt = self.get_pharmvar_allele_info()
            for el in dt:
                dal = dt[el][0]

                gene = dal["geneSymbol"]
                allele = dal["coreAllele"]
                function = dal["function"]
                if( dal["variants"] is not None ):
                    variants = list( filter( lambda x: "rsId" in list(x), dal["variants"] )) 
                    variants = list( map( lambda x: x["rsId"], variants ))

                    if(not gene in dat):
                        dat[gene] = {}
                    if( not allele in dat[gene] ):
                        dat[gene][allele] = {}
                    dat[gene][allele]["function"] = function
                    dat[gene][allele]["variants"] = variants

            json.dump(dat, open(opath, 'w') )
        else:
            dat = json.load( open(opath, 'r') )

        return dat

    def parse_alleleSpecific_mutations_in_cases(self):
        """
Steps:
    - get genes with named alleles in pharmvar
    - get allele functions
    - download and parse allele definition table (clinPgx)
    - parse tcga mutation file of samples and get the mutations and the alleles that they happen, it 
        """
        
        """
        dat = self.get_pharmvar_alleles()
        ids = set([ x["pvId"] for x in dat ])
        genes = set([ x["geneSymbol"] for x in dat ])
        infodat = self.get_pharmvar_allele_info(ids)
        """
        dat = self.process_allele_info()

        # The function below is not necessary because the individual info from the alleles already have the variants
        #defdat = self._download_parse_allele_definition_table(genes)

        for project in tqdm(self.projects):
            indir = os.path.join(self.data_dir, project)
            dcases = self.proc._get_cases_metadata(project)

            datcat = 'simple nucleotide variation'
            odir, fsodir, file_list = self.proc.get_case_files_by_data_category(project, datcat)

            start_alleles_cases = {}

            mapp = self.proc._get_map_case_file(odir)
            mapp_tissue = self.proc._get_map_file_condition(odir)
            for uuid in file_list:
                f = "raw_%s.out" %(uuid)
                path = os.path.join(fsodir, f)
                case_id = mapp[uuid]
                tissue_type = mapp_tissue[uuid]

                if(case_id in dcases):
                    if( not case_id in start_alleles_cases):
                        start_alleles_cases[case_id] = { "mutations": [] }

                    dclin = dcases[case_id]
                    snps = self._get_snprs_annotation(path)

                    gender = dclin["gender"]
                    race = dclin["race"]
                    ethnicity = dclin["ethnicity"]
                    muts = list(snps)
                    for m in muts:
                        gene = m.split('-')[0]
                        if( gene in dat ):
                            refb = snps[m]["reference_allele"]
                            a1b = snps[m]["tumor_seq_allele1"]
                            a2b = snps[m]["tumor_seq_allele2"]
                            rsid = snps[m]["dbsnp_rs"]
                            
                            alleles_mutated = []
                            if( refb!=a1b ):
                                alleles_mutated.append("a1")
                            if( refb!=a2b ):
                                alleles_mutated.append("a2")

                            genotype = f"{refb},{a1b}|{a2b}"
                            start_alleles_cases[case_id]["mutations"].append( { "gene": gene, "gtype": genotype, "rsid": rsid, "alleles_mutated": alleles_mutated } )

            matched_cases = set()
            for case_id in start_alleles_cases:
                start_alleles_cases[case_id]["identified_alleles"] = {}

                genes = set( map( lambda x: x["gene"], start_alleles_cases[case_id]["mutations"] ))

                for g in genes:
                    filtered_case_mutations = list( filter( lambda x: x["gene"] == g, start_alleles_cases[case_id]["mutations"] ))
                    filtered_case_mutations =  list( map( lambda x: x["rsid"], filtered_case_mutations ))
                    #print(case_id, g, filtered_case_mutations)
                    galleles = dat[g]
                    for al in galleles:
                        smuts = galleles[al]["variants"]
                        t = len(smuts)
                        inter = set(filtered_case_mutations).intersection( set(smuts) )
                        n = len(inter)
                        if( n > 0 ):
                            matched_cases.add(case_id)
                            if( not al in start_alleles_cases[case_id]["identified_alleles"] ):
                                start_alleles_cases[case_id]["identified_alleles"][al] = []
                            
                            info = list( filter( lambda x: ( x["gene"] == g and x["rsid"] in inter ), start_alleles_cases[case_id]["mutations"] ))
                            info =  list( map( lambda x: { "rsid": x["rsid"], "genotype": x["gtype"], "alleles_mutated": x["alleles_mutated"]}, info ))

                            allinfo = { "gene": g, "coverage": (n/t), "details": info }
                            start_alleles_cases[case_id]["identified_alleles"][al].append( allinfo )

            print( project, " - cases with some PharmVar allele: ", len(matched_cases) )
            opath = os.path.join( odir, "pgx_alleles_information.json" )
            json.dump( start_alleles_cases, open(opath, 'w') )


    def __check_gene_name(self, chromosome, start, end):
        chromosome = chromosome.lower().replace('chr','')
        start = int(start)
        end = int(end)

        r = []
        keys = list( filter(  lambda x: (self.gtf_map[x]['chrom'] == chromosome) and (self.gtf_map[x]['coord'][0] >= start) and (self.gtf_map[x]['coord'][1] <= end), self.gtf_map ) )
        for x in keys:
            r.append( [ self.gtf_map[x]['gid'], x, self.gtf_map[x]['chrom'], self.gtf_map[x]['coord'][0], self.gtf_map[x]['coord'][1] ] )
        return r

    def get_specific_allele_count(self):
        self.gtf_map = json.load( open('/mnt/yasdata/home/yasmmin/Dropbox/portfolio_2025/gdc_explorer_app_web/gdc_exploration_project/copy_number_variation/out/map_gtf.json') )
        path = '/home/yasmmin/Downloads/gdc_download_20260616_145530.067725/db58c3e8-879c-406d-a37a-ac31d4100b3e/1289b406-41ff-445e-93c4-075d0a2edb25.wgs.ASCAT.copy_number_variation.seg.txt'
        
        d = {}
        df = pd.read_csv(path, sep='\t')
        for i in df.index:
            chrom = df.loc[i, 'Chromosome']
            s = df.loc[i, 'Start']
            e = df.loc[i, 'End']
            r = self.__check_gene_name(chrom, s, e)
            print(r)
            for x in r:
                if( not x[1] in d ):
                    d[ x[1] ] = 0
                d[ x[1] ] += 1

        maxx = 0
        maxg = 0
        for g in d:
            if( d[g] > maxx ):
                maxx = d[g]
                maxg = g
            #print(g, d[g])
        print('max count', d[g], '-', g)

    def explore_found_alleles_project_cases(self):
        # when allele 1 has an empty array, it is the default gene*1 allele, which always has normal function. However individuals that have *1xN or *2xN have increased function, and are rapid or ultra rapid metabolizers.
        # Reference eaning of columns in maf files: https://docs.gdc.cancer.gov/Data/File_Formats/MAF_Format/
        # The columns t_ref_count and t_alt_count mean the depth supporting referene and alternative allele in bam file, so the major and minor allele to integrate with CNV data for copies is basically check who is higher than the other (ref or alt). if Tumor_Seq_Allele1==Tumor_Seq_Allele2 the individual is homozygous
        # The alleles

        idat = self.process_allele_info()
        summary_genes = set()
        variant_pairs = set()

        pjs = self.projects
        pjs = ["TCGA-UCEC"]
        for project in tqdm(self.projects):
            indir = os.path.join(self.data_dir, project)
            dcases = self.proc._get_cases_metadata(project)

            datcat = 'simple nucleotide variation'
            odir, fsodir, file_list = self.proc.get_case_files_by_data_category(project, datcat)
            path = os.path.join( odir, "pgx_alleles_information.json" )
            dat = json.load( open(path, 'r') )
            dc = {}
            for case in dat:
                if( len(dat[ case ]['identified_alleles']) > 0 ):
                    alleles = dat[ case ]['identified_alleles']
                    for allele in alleles:
                        gene = allele.split("*")[0]
                        summary_genes.add(gene)

                        try:
                            function = idat[gene][allele]["function"]
                        except:
                            function = ""
                        if( 'details' in alleles[allele][0] ):
                            for info in alleles[allele][0]["details"]:
                                rsid = info["rsid"]
                                if(rsid!='nan'):
                                    if(not case in dc):
                                        dc[case] = { 'a1': [], 'a2': [] }
                                    gtype = info["genotype"]
                                    for am in info['alleles_mutated']:
                                        dc[case][am].append([gene, allele, function, rsid, gtype])
                                        variant_pairs.add( (gene, rsid) )

            print(project, dc)
            opath = os.path.join( odir, "summary_allele_info.json" )
            json.dump( dc, open( opath, 'w' ) )

        summary_genes = list(summary_genes)
        opath = os.path.join( self.out, "pgx_clingx_summary_genes.json" )
        json.dump( summary_genes, open( opath, 'w' ) )

        variant_pairs = list(variant_pairs)
        opath = os.path.join( self.out, "pgx_clingx_variant_pairs.json" )
        json.dump( variant_pairs, open( opath, 'w' ) )
        
    def _clean_sumary_data(self, g, od):
        dat = {}
        for d in od:
            try:
                genes = list( map( lambda x: x["symbol"], d["location"]["genes"] ))
                drugs = list( map( lambda x: x["name"], d["relatedChemicals"] ))
                diseases = list( map( lambda x: x["name"], d["relatedDiseases"] ))
                for al in d["allelePhenotypes"]:
                    allele = al["allele"]
                    if(allele.startswith("*")): # because it may be in form of CT
                        pheno = al["phenotype"]
                        for g in genes:
                            allele = g+allele
                            if( not allele in dat ):
                                dat[allele] = []
                            dat[allele].append({ "ann": pheno, "drugs": drugs, "diseases": diseases })
            except:
                print(g, d)
        return dat

    def retrieve_clingx_summary_annotations(self):
        path = os.path.join( self.out, "pgx_clingx_summary_genes.json" )
        genes = json.load( open( path, 'r' ) )

        #important metadata: allelePhenotypes.allele, allelePhenotypes.phenotype, relatedChemicals.name, relatedDiseases.name
        dpath = os.path.join(self.out, "summary_information")
        if( not os.path.isdir(dpath) ):
            os.makedirs(dpath)

        dat = {}
        for g in tqdm(genes):
            tmp = {}
            opath = os.path.join( dpath, f"{g.lower()}.json" )
            if( not os.path.exists(opath) ):
                url = f"https://api.clinpgx.org/v1/data/summaryAnnotation?location.genes.symbol={g}"
                r = requests.get( url )
                tmp = r.json()["data"]
                json.dump( tmp, open( opath, 'w' ) )
                sleep(2)
            else:
                tmp = json.load( open( opath, 'r' ) )
            dat[g] = self._clean_sumary_data( g.lower(), tmp )

        return dat

    def _clean_variant_data(self, _id, od):
        dat = []
        for d in od:
            try:
                try:
                    gene = d["geneSymbol"]
                except:
                    gene = _id.split('_')[0].upper()
                rsid = d["location"]["rsid"]

                genotype = d["alleleGenotype"].split(" + ")
                assoc = d["isAssociated"]
                sentence = d["sentence"]
                phenoCats = list( map( lambda x: x["term"], d["phenotypeCategories"] ))
                drugs = list( map( lambda x: x["name"], d["relatedChemicals"] ))
                pops = set()
                if( "populationPhenotypes" in d ):
                    for it in d["populationPhenotypes"]:
                        if( "disease" in it ):
                            pops.add( it["disease"]["name"] )
                pops = list(pops)

                # dat[_id] = { "genotype": genotype, "ann": sentence, "drugs": drugs, "phenoCats": phenoCats, "isAssociated": assoc, "population": pops }
                dat.append({ "genotype": genotype, "ann": sentence, "drugs": drugs, "phenoCats": phenoCats, "isAssociated": assoc, "population": pops })
            except:
                print(_id, d)
        return dat

    def retrieve_clingx_variant_annotations(self):
        path = os.path.join( self.out, "pgx_clingx_variant_pairs.json" )
        pairs = json.load( open( path, 'r' ) )

        # important metadata: alleleGenotype, sentence, isAssociated, phenotypeCategories.term, relatedChemicals.name, populationPhenotypes.disease.name
        # the allele may come as AA + AG (composed) or alone A
        dpath = os.path.join(self.out, "variant_information")
        if( not os.path.isdir(dpath) ):
            os.makedirs(dpath)

        dat = {}
        for p in tqdm(pairs):
            g, v = p
            _id = f"{g.lower()}_{v.lower()}"
            tmp = {}
            opath = os.path.join( dpath, f"{_id}.json" )
            if( not os.path.exists(opath) ):
                url = f"https://api.clinpgx.org/v1/data/variantAnnotation?location.genes.symbol={g}&location.fingerprint={v}"
                r = requests.get( url )
                tmp = r.json()["data"]
                json.dump( tmp, open( opath, 'w' ) )
                sleep(2)
            else:
                tmp = json.load( open( opath, 'r' ) )
            dat[_id] = self._clean_variant_data( _id, tmp )

        return dat

    def get_annotations(self):

        # 1 - get all genes observed in projects and the rsids - ok
        # 2 - with the genes get summary annotations - ok
        # 3 - with the gene-rsid pairs get variant annotations - ok

        # 4 - personalize case report with the specific case genotype of alleles and variants
            # button to modal for summary, and another for variant
            # highlight genotype
            # Cross the annotations regarding the gene & rsid related to the selected mutation?

        sann = self.retrieve_clingx_summary_annotations()
        vann = self.retrieve_clingx_variant_annotations()

        for project in tqdm(self.projects):
            indir = os.path.join(self.data_dir, project)
            dcases = self.proc._get_cases_metadata(project)

            datcat = 'simple nucleotide variation'
            odir, fsodir, file_list = self.proc.get_case_files_by_data_category(project, datcat)
            path = os.path.join( odir, "summary_allele_info.json" )
            if( os.path.exists(path) ):
                dat = json.load( open(path, 'r') )
                # dc[case][am].append([gene, allele, function, rsid, gtype])

                dc = {}
                for case_id in dat:
                    all_alleles = {}

                    info = dat[case_id]
                    vdat = []
                    for am in info:
                        for el in info[am]:
                            gene, allele, function, rsid, gtype = el
                            if( gene in sann ):
                                if( allele in sann[gene] ):
                                    if( not allele in all_alleles ):
                                        all_alleles[allele] = { 'function': function, 'anns': [] }
                                    all_alleles[ allele ]['anns'].append( sann[gene][allele] )

                            _id = f"{ gene.lower() }_{ rsid.lower() }"
                            if( _id in vann ):
                                vd = vann[_id]
                                for it in vd:
                                    vgtype = it["genotype"]
                                    vg1 = gtype.split(",")[-1].replace("|", "")
                                    vg2 = vg1[-1]
                                    if( (vg1 in vgtype) or (vg2 in vgtype) ):
                                        it["gene"] = gene
                                        it["rsid"] = rsid
                                        vdat.append(it)

                    dc[case_id] = { 'summary_clingx': all_alleles, 'variant_clingx': vdat  }
                    #if(len(vdat) > 0):
                    #    print(project, case_id, vdat)

                #print(project, dc)
                opath = os.path.join( odir, "clingx_filtered_annotation_info_for_cases.json" )
                json.dump( dc, open( opath, 'w' ) )

                # ---> main div : analysis_clinpgx_ann_ai3

                # cases_summary_table_ai3_clinpgx
                #    table clinpgx_summary_datatab - cols: case, allele, function, diseases (badge dark blue), drugs (badge light blue), sentence (tag p small text size)

                # cases_variant_table_ai3_clinpgx
                #   table clinpgx_variant_datatab - cols: case, variant (rsid in gene hugo_symbol), genotype, phenoCategory, drugs, sentence

    def integrate_cnv_n_allele_copies(self):
        # Only coad has cnv data stored, and the cnv files contain the major copy allele copy count and the same for minor, if minor equals to 0, it refers to a loss of heterozygosity
        # study case for case_id e703af88-c05f-4d7b-b225-d640e4d1a2a5
        snv_path = "/var/www/html/gdc_explorer_app/data_processed/TCGA-COAD/simple-nucleotide-variation/files/raw_39495aca-acc5-417a-8679-22e9246e77c0.out"
        cnv_path = "/var/www/html/gdc_explorer_app/data_processed/TCGA-COAD/copy-number-variation/files/raw_8baa452a-356d-4765-8ab2-c5c320ed7f78.out"

        cdf = pd.read_csv(cnv_path, sep='\t', comment="#")
        cdf = cdf.dropna()
        genes = list( map( lambda x: x.split('.')[0], cdf[ (cdf["min_copy_number"] != 1) & (cdf["max_copy_number"] != 1 )  ]["gene_name"].values ))
        #print(genes)
        sdf = pd.read_csv(snv_path, sep='\t', comment="#")
        f = sdf[ (sdf["Hugo_Symbol"].isin(genes)) & ~(sdf["dbSNP_RS"].isna()) ][ ['Reference_Allele', 'Tumor_Seq_Allele1', 'Tumor_Seq_Allele2', 'HGVSc', 'Hugo_Symbol', 'SWISSPROT', 'Variant_Classification', 'Consequence', 'IMPACT', 'VARIANT_CLASS', 'dbSNP_RS', 'SIFT', 'PolyPhen', "t_ref_count", "t_alt_count"] ]
        print(f)

        f2 = sdf[ sdf["Tumor_Seq_Allele1"] != sdf["Reference_Allele"] ] # zero count - checking if there is a possibility of al1 and al2 being mutated, and maybe for a distinct base

        f2 = sdf[ ['Reference_Allele', 'Tumor_Seq_Allele1', 'Tumor_Seq_Allele2', 'HGVSc', 'Hugo_Symbol', 'SWISSPROT', 'Variant_Classification', 'Consequence', 'IMPACT', 'VARIANT_CLASS', 'dbSNP_RS', 'SIFT', 'PolyPhen', "t_ref_count", "t_alt_count"] ][ sdf["t_ref_count"] < sdf["t_alt_count"] ]
        """
        There are two entries, one it is for smad4 gene, but it has no dbsnp id associated, but for the same position and gene here is a entry in dbsnp (rs2144452610) but changing from G to A or C, and not for T.
        The read count is very low for both columns (tref=8, and talt=15). Since it is one of the major tumor supressor genes, it is frequently mutated and inactivated (PMC8054659) so this could explain the underexpression
        """
        
        print(f2)

    def run(self):
        #self.merge_datasources()
        #self.generate_app_pgx_compiled_data()
        #self._prepare_drugs_list()
        # treat drugnames because they were annotated composed (+)
        #self.get_smiles_from_graphql_pdb()

        #self.clean_drug_names()
        #self.get_clinpgx_drugLabels()

        #self.get_pharmvar_alleles()
        
        #self.parse_alleleSpecific_mutations_in_cases()
        #self.process_allele_info()
        
        #self.explore_found_alleles_project_cases()
        #self.retrieve_clingx_summary_annotations()
        #self.retrieve_clingx_variant_annotations()

        self.get_annotations()

        #self.integrate_cnv_n_allele_copies()

if( __name__ == "__main__" ):
    out = '/mnt/yasdata/home/yasmmin/Dropbox/portfolio_2025/gdc_explorer_app_web/gdc_exploration_project/bridge_pharcoGenomics/out'
    
    o = HandleBridgePharcogenomicsAnalysis( out)
    o.run()