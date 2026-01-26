import os
import json
import gzip

import numpy as np
import pandas as pd
import gseapy as gp

from tqdm import tqdm
from scipy import stats

from process_data import DataWrangler

class HandleMethylationAnalysis:
    def __init__(self, dir_gdc_cpg_map, fout):
        self.cpg_map = dir_gdc_cpg_map

        self.out = fout
        if( not os.path.exists(self.out) ):
            os.makedirs( fout )

        self.proc = DataWrangler()

    def _make_cpg_mapping(self):
        cpgs = {}
        opath = os.path.join( self.out, 'all_map.json' )
        if( not os.path.isfile(opath) ):
            for f in os.listdir( self.cpg_map ):
                if( f.endswith('.gz') ):
                    path = os.path.join( self.cpg_map, f )
                    with gzip.open( path, 'rb' ) as g:
                        i = 0
                        for line in g:
                            if( i>0 ):
                                l = line.decode('utf-8').split('\t')

                                cpg = l[4]
                                chrom = l[0]
                                start = l[1]
                                end = l[2]
                                gene = l[5]
                                _id = "%s_%s_%s_%s" %( chrom, gene.split(';')[0], start, end )

                                cpgs[ cpg ] = { 'id': _id, 'chromosome': chrom, 'start': start, 'end': end, 'gene': gene }
                            i += 1
            json.dump( cpgs, open(opath, 'w') )
        else:
            cpgs = json.load( open(opath, 'r') )

        return cpgs

    def _get_covariables_table_glm(self, odir, map_uuid_treat, cmeta):
        mp = self.proc._get_map_case_file(odir)
        cols_stratification = ['race','gender', 'ethnicity']
        map_covs = { 'condition': { } }
        for c in cols_stratification:
            map_covs[c] = {}

        body = []
        for uid in map_uuid_treat:
            covs = []
            condition = map_uuid_treat[uid]
            if( not condition in map_covs['condition'] ):
                map_covs['condition'][condition] = len( map_covs['condition'] ) + 1
            covs.append( map_covs['condition'][condition] )

            case_id = mp[uid]
            dat = cmeta[case_id]
            for c in cols_stratification:
                v = dat[c]
                if( not v in map_covs[c] ):
                    map_covs[c][v] = len( map_covs[c] ) + 1
                covs.append( map_covs[c][v] )

            el = [ "%s_%s" %(condition, uid) ] + covs
            body.append( el )

        header = [ ['sampleID'] + cols_stratification ]
        lines = header + body
        tpath = os.path.join(odir, 'table_covariables.tsv' )
        ls = list( map( lambda x: ','.join( [ str(y) for y in x ] ), lines ))
        f = open(tpath, "w")
        f.write("\n".join(ls) + "\n")
        f.close()

        return tpath

    def _get_table_betaValues(self, odir, map_uuid_treat, metdat, filtered_cpgs):
        header = ['probeID']
        for uid in map_uuid_treat:
            condition = map_uuid_treat[uid]
            header.append( "%s_%s" %(condition, uid) )

        body = []
        for cpg in filtered_cpgs:
            values = [cpg]
            for uid in map_uuid_treat:
                try:
                    v = metdat[cpg][uid]
                except:
                    v = 0
                values.append( v )
            body.append(values)

        lines = header + body
        tpath = os.path.join(odir, 'table_databeta.tsv' )
        ls = list( map( lambda x: ','.join( [ str(y) for y in x ] ), lines ))
        f = open(tpath, "w")
        f.write("\n".join(ls) + "\n")
        f.close()

        return tpath

    def _prioritize_genes_cpgs(self, map_cpg, metdat):
        # check genes related to colon mucinous adenocarcinoma
        incommon = set(map_cpg).intersection( set(metdat) )
        genes = set( )
        genes.update( [ map_cpg[x]['gene'].split(';') for x in incommon ] )
        print('total genes', len(genes) )

        # Choosing by those genes that have the highest variance across the samples
        values = {}
        var = {}
        for k in incommon:
            values[k] = [ x for x in metdat[k].values() ]
            var[k] = np.var( values[k] )
        
        revar = dict( sorted( var.items(), key=lambda item: item[1], reverse=True) )
        top = list(revar)[:200]
        tgenes = set( [ map_cpg[x]['gene'].split(';')[0] for x in top ] )

        # selection by those that are associated to any cancer pathway
        ags = set( [ map_cpg[x]['gene'].split(';')[0] for x in list(revar) ] )
        subset = 'c6'
        msig = gp.Msigdb()
        gmt = msig.get_gmt(category='%s.all' %(subset), dbver="2025.1.Hs")

        enr = gp.enrichr(gene_list = list(ags), gene_sets = gmt, organism='human', outdir=None )
        df = enr.res2d
        fgs = set()
        for gs in df['Genes']:
            fgs.update( gs.split(';') )
        filtered_cpgs = list( filter( lambda x: map_cpg[x]['gene'] in fgs, incommon )) # 161528

        return filtered_cpgs

    def _get_exp(self):
        project = 'TCGA-COAD'
        datcat = 'dna methylation'
        odir, fsodir, file_list = self.proc.get_case_files_by_data_category(project, datcat)
        ok_cases, map_uuid_treat = self.proc._select_cases_with_normalAndTumor(odir)
        cpgs = self._make_cpg_mapping()
        cmeta = self.proc.get_cases_info_by_project(project)
        grpcov_path = self._get_covariables_table_glm(odir, map_uuid_treat, cmeta)

        path = os.path.join(odir, 'parsable_cpgs_info.json')
        metdat = json.load( open(path, 'r') )

        filtered_cpgs = self._prioritize_genes_cpgs( cpgs, metdat)
        data_path = self._get_table_betaValues( odir, map_uuid_treat, metdat, filtered_cpgs)

        # cpgtools call
        opath = os.path.join( odir, "GLM_met_%s" %( project.lower() ) )
        os.system( "dmc_glm.py  -i %s -g %s -o %s" %(data_path, grpcov_path, project.lower() ) )
        # dmc_glm.py  -i test_05_TwoGroup.tsv.gz -g test_05_TwoGroup.grp.csv -o GLM_2G

    def run(self):
        self._make_cpg_mapping()

if( __name__ == "__main__" ):
    cpg_map = '/mnt/yasdata/home/yasmmin/Dropbox/iarc_job/gdc_exploration_project/methylation/gdc_cpg_mapping/'
    out = '/mnt/yasdata/home/yasmmin/Dropbox/iarc_job/gdc_exploration_project/methylation/out'
    
    cpg_map = '/aloy/home/ymartins/cigs_exercise/gdc_exploration_project/methylation/gdc_cpg_mapping/'
    out = '/aloy/home/ymartins/cigs_exercise/gdc_exploration_project/methylation/gdc_cpg_mapping/out'
    o = HandleMethylationAnalysis(cpg_map, out)
    o.run()