import os
import json
import gzip

import numpy as np
import pandas as pd
import gseapy as gp

from tqdm import tqdm
from scipy import stats

from process_data import DataWrangler

# https://tcpa.drbioright.org/rppa500/differential_analysis.html

class HandleProteomicsAnalysis:
    def __init__(self, fout):
        self.compiled_tcga = '../external_db/compiled_proteomic_tcga.tsv'

        self.out = fout
        if( not os.path.exists(self.out) ):
            os.makedirs( fout )

        self.proc = DataWrangler()

    def _get_exp(self):
        project = 'TCGA-COAD'
        datcat = 'proteome profiling'
        odir, fsodir, file_list = self.proc.get_case_files_by_data_category(project, datcat)
        print( "Number of files", len(file_list) )
        
        # Retrieving raw files
        for uuid in file_list:
            self.proc._get_file_by_uuid( fsodir, uuid)

    def _load_mapping_submitId_uuid(self, odir):
        path = os.path.join(odir, "files_metadata.tsv")
        df = pd.read_csv(path, sep='\t')
        dc = dict( zip( df['file_id'].values, df['cases.0.submitter_id'].values ))

        return dc

    def _check_match_compiled_vs_rawRppa(self):
        project = 'TCGA-COAD'
        datcat = 'proteome profiling'
        odir, fsodir, file_list = self.proc.get_case_files_by_data_category(project, datcat)
        
        mp = self._load_mapping_submitId_uuid(odir)
        df = pd.read_csv( self.compiled_tcga, sep='\t')
        genes = list(df.columns)[1:]
        gene_in_samples = set()
        stats = {}
        ratios = []
        dat = {}
        for uuid in tqdm(file_list):
            f = "raw_%s.out" %(uuid)
            submitId = mp[uuid]
            ref = dict( zip( genes, df[ df.sample_id.str.contains(submitId) ][genes].values[0] ))
            sampleId = df[ df.sample_id.str.contains(submitId) ]['sample_id'].values[0]

            path = os.path.join(fsodir, f)
            idf = pd.read_csv(path, sep='\t')
            idf = idf[ idf.peptide_target.isin(genes) ]
            igenes = idf.peptide_target.values
            gene_in_samples.update(igenes.tolist())

            alt = dict( zip( igenes, idf.protein_expression.values ) )
            cnt = 0
            for g in igenes:
                if( ref[g] == alt[g] ):
                    cnt += 1

            ttype = 'primary'
            if(cnt == 0):
                ttype = 'metastatic'
            dat[uuid] = { "sample_id": sampleId, "tumor_type": ttype, "expression": alt }

            # from those that are in the main table, the ratio measures how many matched with the value in the table
            ratio = cnt/len(igenes)
            stats[uuid] = [ submitId, sampleId, cnt, ratio, len(genes), len(igenes)  ]
            ratios.append(ratio)
        
        opath = os.path.join(odir, "matching_report.tsv")
        f = open(opath, 'w')
        f.write("uuid\tsubmitter_id\tsample_id\tcount_matches\tratio_matches\ttotal_genes_compiled\ttotal_genes_instance\n")
        for k in stats:
            el = [k] + stats[k]
            line = '\t'.join( [ str(x) for x in el ] )
            f.write(line+'\n')
        f.close()

        ratios = list( filter( lambda x: x > 0, ratios ))
        print('Min match ratio:', min( ratios ) ) # 0.88
        print('Avg match ratio:', sum( ratios )/len( ratios ) ) # 0.930
        print('Max match ratio:', max( ratios ) ) # 0.933

        opath = os.path.join(odir, "all_genes_in_project.json")
        obj = { "genes": list(gene_in_samples) }
        json.dump( obj, open(opath, 'w') )
        
        opath = os.path.join(odir, "prot_exp_project.json")
        obj = { "genes": list(gene_in_samples) }
        json.dump( dat, open(opath, 'w') )
        
        '''
        There were three cases in which none of the values matched that belongs to metastatic tumrs instead of primary tumor
        '''

    def run(self):
        #self._get_exp()
        self._check_match_compiled_vs_rawRppa()

if( __name__ == "__main__" ):
    out = '/mnt/yasdata/home/yasmmin/Dropbox/portfolio_2025/gdc_explorer_app_web/gdc_exploration_project/proteomic/out'
    
    o = HandleProteomicsAnalysis( out)
    o.run()