import os
import json
import gzip

import numpy as np
import pandas as pd
import gseapy as gp

from tqdm import tqdm
from scipy import stats

from process_data import DataWrangler

# cpgtools reference: https://cpgtools.readthedocs.io/en/latest/demo/dmc_glm.html
'''
Check the coefficients importance of each covariable
'''

class HandleCNAAnalysis:
    def __init__(self,  fout):

        self.out = fout
        if( not os.path.exists(self.out) ):
            os.makedirs( fout )

        self.proc = DataWrangler()

    def _get_exp(self):
        project = 'TCGA-COAD'
        datcat = 'copy number variation'
        odir, fsodir, file_list = self.proc.get_case_files_by_data_category(project, datcat)
        ok_cases, map_uuid_treat = self.proc._select_cases_with_normalAndTumor(odir)
        
    def run(self):
        self._get_exp()

if( __name__ == "__main__" ):
    out = '/mnt/yasdata/home/yasmmin/Dropbox/portfolio_2025/gdc_explorer_app_web/gdc_exploration_project/copy_number_variation/out'
    
    o = HandleMethylationAnalysis( out)
    o.run()