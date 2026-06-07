import os
import json
import gzip
import shutil

import numpy as np
import pandas as pd
import gseapy as gp
import statistics as st

from sklearn.feature_selection import RFECV
from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import AdaBoostClassifier
from sklearn.model_selection import cross_val_score

from tqdm import tqdm
from scipy import stats
from scipy.stats import pearsonr
from subprocess import run, PIPE

from biomart import BiomartServer

from process_data import DataWrangler

# cpgtools reference: https://cpgtools.readthedocs.io/en/latest/demo/dmc_glm.html
'''
Check the coefficients importance of each covariable
'''

class HandleCNAAnalysis:
    def __init__(self, fout, gtf = None ):
        self.ctypes = ['del', 'loss', 'neutral', 'gain', 'amp']

        self.out = fout
        if( not os.path.exists(self.out) ):
            os.makedirs( fout )

        self.gtf_map = None
        if( (gtf is not None) and (os.path.exists(gtf)) ):
            self.gtf_map = self.process_gtf(gtf)

        self.out_processed = os.path.join(fout, 'processed')
        if( not os.path.exists(self.out_processed) ):
            os.makedirs( self.out_processed )

        self.proc = DataWrangler()

        #server = BiomartServer("http://www.ensembl.org/biomart")
        #self.dataset = server.datasets['hsapiens_gene_ensembl']

    
    def process_gtf(self, path):
        mp = {}
        opath = os.path.join( self.out, 'map_gtf.json')
        if( not os.path.exists(opath) ):
            f = open(path, 'r')
            for line in f:
                if( not line.startswith('#') ):
                    l = line.replace('\n','').split('\t')
                    if( l[2] == 'CDS' ):
                        chrom = l[0]
                        start = int(l[3])
                        end = int(l[4])
                        info = l[-1]

                        gid = None
                        gname = None
                        gversion = '.'
                        for g in info.split('; '):
                            if( g.startswith('gene_id') ):
                                gid = g.replace('gene_id ','').replace('"','')
                            if( g.startswith('gene_name') ):
                                gname = g.replace('gene_name ','').replace('"','')
                            if( g.startswith('gene_version') ):
                                gversion += g.replace('gene_version ','').replace('"','')
                                gid += gversion

                        if( (gid != None) and (gname != 'None') ):
                            if( not gname in mp ):
                                mp[gname] = { 'gid': gid, 'chrom': chrom, 'coord': [start, end] }
                            else:
                                if( start < mp[gname]['coord'][0] ):
                                    mp[gname]['coord'][0] = start
                                if( end > mp[gname]['coord'][1] ):
                                    mp[gname]['coord'][1] = end
            f.close()

            json.dump( mp, open(opath, 'w') )
        else:
            mp = json.load( open(opath, 'r') )

        return mp

    def __check_gene_name(self, chromosome, start, end):
        r = []
        
        chromosome = chromosome.lower().replace('chr','')
        start = int(start)
        end = int(end)

        if( self.gtf_map is not None ):
            keys = list( filter(  lambda x: (self.gtf_map[x]['chrom'] == chromosome) and (self.gtf_map[x]['coord'][0] >= start) and (self.gtf_map[x]['coord'][1] <= end), self.gtf_map ) )
            for x in keys:
                r.append( [ self.gtf_map[x]['gid'], x, self.gtf_map[x]['chrom'], self.gtf_map[x]['coord'][0], self.gtf_map[x]['coord'][1] ] )
        else:
            response = self.dataset.search({
                'filters': {
                    'chromosome_name': chromosome,
                    'start': start,
                    'end': end
                },
                'attributes': [
                    'ensembl_gene_id',
                    'external_gene_name',
                    'chromosome_name',
                    'start_position',
                    'end_position'
                ]
            })

            for line in response.iter_lines():
                r.append( line.decode('utf-8').split('\t') )

        return r

    def _get_absolute_copy_number_from_segmentMean(self, path):
        # https://cnvkit.readthedocs.io/en/stable/pipeline.html#call

        df = pd.read_csv(path, sep='\t')
        odf = df[ ['Chromosome','Start','End','Segment_Mean'] ]
        odf['gene'] = [ 'gene'+str(i) for i in range(len(df)) ]
        odf = odf[ ['Chromosome','Start','End','gene','Segment_Mean'] ]
        odf.columns = ['chromosome','start','end', 'gene', 'log2']
        op = 'test_cnv.cns'
        odf.to_csv(op, sep='\t', index=None)

        os.system('cnvkit.py call %s -o tmp.cns' %(op) )
        ip = 'tmp.cns'
        df = pd.read_csv(ip, sep='\t')
        acnvs = df.cn.values.tolist()

        os.remove(op)
        os.remove(ip)

        return acnvs

    def standardize_gdc_cnv_file(self, cond, path):
        """
        Function to standardize the copy number files. 
        In case the pool of files of the samples are mixed, with some having segment_mean and others having the absolute copy number
        """

        max_cols = ['major_copy_number', 'max_copy_number']
        min_cols = ['minor_copy_number', 'min_copy_number']
        
        header = ['gene_id', 'gene_name', 'chromosome', 'start', 'end', 'copy_number', 'min_copy_number', 'max_copy_number']

        fname = "%s_%s" %( cond, path.split( os.path.sep )[-1] )
        opath = os.path.join(self.out_processed, fname)

        df = pd.read_csv(path, sep='\t')
        cols = [ x.lower() for x in df.columns ]
        df.columns = cols

        if( (not os.path.exists(opath)) ):
            if(  ('gene_name' not in cols) ):
                if( 'copy_number' not in cols ):
                    # Corresponds to the datatypes: Masked Copy Number Segment and 
                    acnvs = self._get_absolute_copy_number_from_segmentMean(path)
                else:
                    acnvs = df.copy_number.values.tolist()
                
                lstr = set(  )
                j = 0
                for i in df.index:
                    chrom = df.loc[i, 'chromosome']
                    start = df.loc[i, 'start']
                    end = df.loc[i, 'end']

                    min_copy = -1
                    for mc in min_cols:
                        if(mc in df.columns):
                            min_copy = df.loc[i, mc]

                    max_copy = -1
                    for mc in max_cols:
                        if(mc in df.columns):
                            max_copy = df.loc[i, mc]

                    genes = self.__check_gene_name(chrom, start, end)
                    genes = list( filter( lambda x: x[1] != '', genes ))

                    acnv = acnvs[j]

                    for ginfo in genes:
                        x = ginfo + [ acnv, min_copy, max_copy ]
                        x = '\t'.join( [ str(y) for y in x ] )
                        if( x not in lstr ):
                            lstr.add(x)

                    j+=1
                
                #lines = list( map( lambda x: '\t'.join( [ str(y) for y in x ] ), lines ))
                x = '\t'.join( [ str(y) for y in header ] )
                lines = [ x ] + list(lstr)
                f = open( opath, "w")
                f.write("\n".join(lstr) + "\n")
                f.close()
            else:
                mic = 'min_copy_number'
                for mc in min_cols:
                    if(mc in df.columns):
                        mic = mc
                
                mac = 'max_copy_number'
                for mc in max_cols:
                    if(mc in df.columns):
                        mac = mc

                df = df[ [ 'gene_id', 'gene_name', 'chromosome', 'start', 'end', 'copy_number', mic, mac ] ]
                df.columns = header
                df.to_csv( opath, sep='\t', index=None )
    
    def _check_absolute_copy_number_value(self, fsodir, map_uuid_treat):
        acnv = { 'Tumor': 0, 'Normal': 0 }
        for uuid in tqdm(map_uuid_treat):
            cond = map_uuid_treat[uuid]
            path = os.path.join( fsodir, "raw_%s.out" %(uuid) )
            self.standardize_gdc_cnv_file(cond, path)

            df = pd.read_csv(path, sep='\t')
            cols = [ x.lower() for x in df.columns ]
            df.columns = cols
            if( 'copy_number' in cols ):
                acnv[cond] += 1
        print(acnv)

    def extract_common_features(self, odir, fsodir, map_case, map_uuid_treat, meta):
        # Thresholds: https://forum.depmap.org/t/defining-deep-deletions-and-amplifications/710/4

        featdir = os.path.join(odir, "features", "basic")
        if( not os.path.exists(featdir) ):
            os.makedirs(featdir)

        all_feat = set()
        fd = {  } # tumor is 1
        opath = os.path.join(featdir, "features_cna.json")
        if( not os.path.exists(opath) ):
            for uuid in tqdm(map_uuid_treat):
                cond = map_uuid_treat[uuid]
                fname = "%s_raw_%s.out" %( cond, uuid )
                path = os.path.join(self.out_processed, fname)

                try:
                    case = map_case[uuid]
                    _id = "%s_%s" %(case, uuid)
                    
                    fd[_id] = { 'y': 0, 'y_vital': 0 }
                    if(cond == 'Tumor'):
                        fd[_id]['y'] = 1

                    if( 'vital_status' in meta[case] ):
                        vital_status = meta[case]['vital_status']
                        if(vital_status == 'Dead'):
                            fd[_id]['y_vital'] = 1

                    df = pd.read_csv(path, sep='\t')
                    df = df[ ~df.copy_number.isna() ]
                    for i in df.index:
                        gene = df.loc[i, 'gene_name']
                        chrom = df.loc[i, 'chromosome']
                        cnv = int( df.loc[i, 'copy_number'] )
                        class_cnv = 'neutral'
                        if( cnv >= 2.64 and cnv <= 3.36 ):
                            class_cnv = 'gain'
                        elif( cnv > 3.36 ):
                            class_cnv = 'amp'
                        elif( cnv >= 0.87 and cnv <= 1.32 ):
                            class_cnv = 'loss'
                        elif( cnv < 0.87 ):
                            class_cnv = 'del'
                        feat = "%s_%s_%s" %(class_cnv, gene, chrom)
                        all_feat.add(feat)
                        fd[_id][feat] = cnv
                except:
                    print(path)

                json.dump( fd, open( opath,'w') )
        else:
            fd = json.load( open(opath, 'r') )
            for k in fd:
                all_feat.update( list( filter( lambda x: (x.split('_')[0] in self.ctypes) , fd[k].keys() ) ) )

        all_feat = list( filter( lambda x: (x.split('_')[0] in self.ctypes) , all_feat ) )
        print( '# Features:', len(all_feat) )
        header = ['sample', 'y', 'y_vital'] + list(all_feat)
        #lines = [ header ]
        x = '\t'.join( [ str(y) for y in header ] )
        opath = os.path.join(featdir, "feat_table.tsv")
        f = open( opath, "w")
        f.write( x + "\n")
        f.close()
        for k in fd:
            fvs = [ fd[k]['y'], fd[k]['y_vital'] ]
            for ft in all_feat:
                try:
                    v = fd[k][ft]
                except:
                    v = -1
                fvs.append(v)
            el = [k] + fvs
            x = '\t'.join( [ str(y) for y in el ] )
            with open( opath, "a") as f:
                f.write( x + "\n")
            #lines.append(el)

        #lines = list( map( lambda x: '\t'.join( [ str(y) for y in x ] ), lines ))
        #f = open( opath, "w")
        #f.write("\n".join(lstr) + "\n")
        #f.close()

    def _reconstruct_feat_cna_json(self, featdir):
        opath = os.path.join(featdir, "features_cna.json")
        if( not os.path.isfile(opath) or os.path.getsize(opath) < 10 ):
            fd = {}
            path = os.path.join(featdir, "feat_table.tsv")
            i = 0
            f = open(path, 'r')
            for line in f:
                l = line.replace('\n','').split('\t')
                if( i == 0 ):
                    header = l[1:]

                else:
                    dc = {}
                    j = 0
                    tmp = l[1:]
                    for v in tmp:
                        if( v != '-1' ):
                            h = header[j]
                            try:
                                dc[h] = float(v)
                            except:
                                dc[h] = v
                        j += 1
                    fd[ l[0] ] = dc
                i+=1
            f.close()

            json.dump( fd, open(opath, 'w') )

    def analysis_cna_events(self, featdir):
        opath = os.path.join(featdir, "features_cna.json")
        d = json.load( open( opath,'r') )

        # How popular is each feature among the samples?
        m = {}
        for k in d:
            cond = 'Tumor'
            if( d[k]['y'] == 0 ):
                cond = 'Normal'

            if( not cond in m ):
                m[cond]={}
            for f in d[k]:
                if( not f.startswith('y') ):
                    if( not f in m[cond] ):
                        m[cond][f] = 0
                    m[cond][f] += 1

        for cond in m:
            tmp = m[cond]
            m[cond] = dict( sorted( tmp.items(), key=lambda item: item[1], reverse=True ) )

        # getting the most popular from each type
        ots = list( filter( lambda x: (x.split('_')[0] in self.ctypes), list(m['Tumor']) )) # 148925
        ons = list( filter( lambda x: (x.split('_')[0] in self.ctypes), list(m['Normal']) )) # 137091

        print( 'Exclusive in Tumor:', len( set(ots) - set(ons) ) ) # 13978
        print( 'Exclusive in Normal:', len( set(ons) - set(ots) ) ) # 2144
        print( 'Both:', len( set(ots).intersection( set(ons) ) ) ) # 134947

        for c in self.ctypes:
            print( c, ' in Tumor:', len( list( filter( lambda x: x.startswith(c), ots )) ) ) # 78407
            print( c, ' in Normal:', len( list( filter( lambda x: x.startswith(c), ons )) ) ) # 78397

        
        # Plot frequencies?

    def _prepare_input_significance_tes(self, odir):
        featdir = os.path.join(odir, "features", "basic")
        inpath = os.path.join(featdir, "features_cna.json")
        d = json.load( open(inpath, 'r') )

        # Prepare input
        opath = os.path.join(featdir, "input_significance_test.tsv")
        if(not os.path.exists(opath)):
            fg = open( opath, 'w')
            fg.write('feature\tcondition\tcase\tuuid\tcnv_value\n')
            for k in d:
                case, uuid = k.split('_')
                cond = 'Tumor'
                if( d[k]['y'] == 0 ):
                    cond = 'Normal'

                for f in d[k]:
                    v = d[k][f]
                    if( (not f.startswith('y')) and (not f.startswith('sample')) and (not f.startswith('neutral')) ):
                        fg.write( '%s\t%s\t%s\t%s\t%i\n' %(f, cond, case, uuid, v) )

            fg.close()

        ps = {}
        df = pd.read_csv(opath, sep='\t')
        for i in df.index:
            f = df.loc[i, 'feature']
            cond = df.loc[i, 'condition']
            case = df.loc[i, 'case']
            v = df.loc[i, 'cnv_value']
            if( not f in ps ):
                ps[f]={}
            if( not cond in ps[f] ):
                ps[f][cond]={}
            if( not case in ps[f][cond] ):
                ps[f][cond][case] = []
            ps[f][cond][case].append(v)
        
        return df, ps

    def test_significance_cnv_values(self, odir):
        # Prepare for paired ttest - normal vs tumor
        # group the values for del with min, and the amplification with max.
        # enrich with the databases
        
        featdir = os.path.join(odir, "features", "basic")
        ups = ['amp', 'gain']
        downs = ['del', 'loss']

        df, dat = self._prepare_input_significance_tes(odir)
        
        opath = os.path.join(featdir, "result_significance_test.tsv")
        fg = open(opath, 'w')
        fg.write('feature\tnumber_normal_cases\tnumber_tumor_cases\tstats\tp_value\n')
        for f in tqdm(dat):
            ftype = f.split('_')[0]
            if( ('Tumor' in dat[f]) and ('Normal' in dat[f]) ):
                if( ftype in ups ):
                    gn = list(map( lambda x: max(x), dat[f]['Normal'].values() ))
                    gt = list(map( lambda x: max(x), dat[f]['Tumor'].values() ))

                if( ftype in downs ):
                    gn = list(map( lambda x: min(x), dat[f]['Normal'].values() ))
                    gt = list(map( lambda x: min(x), dat[f]['Tumor'].values() ))
                
                if( (len(gn) > 2) and (len(gt) > 2) ):
                    t_stat, p_value = stats.ttest_ind(gn, gt, equal_var=False)
                    fg.write("%s\t%i\t%i\t%.6f\t%.6f\n" %( f, len(gn), len(gt), t_stat, p_value))

        """
        ups = ['amp', 'gain']
        grp_amp = df[ (df['feature'].str.startswith('amp')) | (df['feature'].str.startswith('gain')) ].groupby(['feature','condition','case']).max().reset_index()
        for feat in tqdm(grp_amp.feature.unique()):
            gn = grp_amp[ (grp_amp['condition'] == 'Normal') & (grp_amp['feature'] == feat) ]['cnv_value'].values
            gt = grp_amp[ (grp_amp['condition'] == 'Tumor') & (grp_amp['feature'] == feat) ]['cnv_value'].values
            t_stat, p_value = stats.ttest_ind(gn, gt, equal_var=False)
            fg.write("%s\t%i\t%i\t%.6f\t%.6f\n" %(feat, len(gn), len(gt), t_stat, p_value))

        
        grp_del = df[ (df['feature'].str.startswith('del')) | (df['feature'].str.startswith('loss')) ].groupby(['feature','condition','case']).min().reset_index()
        for feat in tqdm(grp_del.feature.unique()):
            gn = grp_del[ (grp_amp['condition'] == 'Normal') & (grp_amp['feature'] == feat) ]['cnv_value'].values
            gt = grp_del[ (grp_del['condition'] == 'Tumor') & (grp_amp['feature'] == feat) ]['cnv_value'].values
            t_stat, p_value = stats.ttest_ind(gn, gt, equal_var=False)
            fg.write("%s\t%i\t%i\t%.6f\t%.6f\n" %(feat, len(gn), len(gt), t_stat, p_value))
        """

        fg.close()

    def _prepare_input_cnsistent(self, fsodir, map_uuid_treat, map_case, featdir):
        opath = os.path.join(featdir, "input_cnsistent.tsv")
        
        if( not os.path.exists(opath) ):
            # prepare input
            header = "sample_id\tchromosome\tstart_pos\tend_pos\tcn\tminor_cn\tmajor_cn\n"
            fg = open(opath, 'w')
            fg.write(header)
            for uuid in tqdm(map_uuid_treat):
                cond = map_uuid_treat[uuid]
                path = os.path.join( fsodir, "raw_%s.out" %(uuid) )
                case = map_case[uuid]
                _id = "%s_%s_%s" %(cond, case, uuid)
                
                df = pd.read_csv(path, sep='\t')
                cols = [ x.lower() for x in df.columns ]
                df.columns = cols
                if( 'copy_number' in cols ):
                    max_cols = ['major_copy_number', 'max_copy_number']
                    min_cols = ['minor_copy_number', 'min_copy_number']
                    
                    for i in df.index:
                        chrom = df.loc[i, 'chromosome']
                        start = df.loc[i, 'start']
                        end = df.loc[i, 'end']
                        cnv = df.loc[i, 'copy_number']

                        min_copy = ""
                        for mc in min_cols:
                            if(mc in df.columns):
                                min_copy = df.loc[i, mc]

                        max_copy = ""
                        for mc in max_cols:
                            if(mc in df.columns):
                                max_copy = df.loc[i, mc]
                                
                        line = f"{_id}\t{chrom}\t{start}\t{end}\t{cnv}\t{min_copy}\t{max_copy}\n"
                        fg.write(line)
            fg.close()

        #df = pd.read_csv(opath, sep='\t')

        #return df
    
    def _partitionate_input_cnsistent(self, featdir):
        featdirp = os.path.join(featdir, "parts_input_cnsistent")
        if( not os.path.exists(featdirp) ):
            os.makedirs(featdirp)

            # cnsistent acceps or the unique copy number value or the pair of minor & major copy number value, not both
            chunk = 20
            parts = {}
            path = os.path.join(featdir, "input_cnsistent.tsv")
            idx = 1
            ant = ""
            i = 0
            f = open(path, 'r')
            for line in f:
                l = line.replace('\n','').split('\t')
                line = '\t'.join(l[:4]+l[-2:])

                if(i==0):
                    header = line
                else:
                    sample = line.split('\t')[0].split('_')[1]
                    if( sample != ant ):
                        if( (len(parts) % chunk == 0) and (len(parts) > 0) ) :
                            lines = list( map( lambda x: '\n'.join(x), parts.values() ))
                            lines = '\n'.join(lines)
                            opath = os.path.join(featdirp, "input_cnsistent_part%i.tsv" %(idx) )
                            g = open(opath, 'w')
                            g.write(header+"\n")
                            g.write(lines+"\n")
                            g.close()
                            
                            parts = {}
                            idx += 1
                        ant = sample


                    if(sample not in parts):
                        parts[sample] = []
                    parts[sample].append(line)
                i += 1

            if( len(parts) > 0):
                lines = list( map( lambda x: '\n'.join(x), parts.values() ))
                lines = '\n'.join(lines)
                path = os.path.join(featdir, "input_cnsistent_part%i.tsv" %(idx) )
                g = open(opath, 'w')
                g.write(header+"\n")
                g.write(lines+"\n")
                g.close()
            f.close()
    
    def _partitionate_input_cnsistent_clean_missing_neutral(self, featdir):
        featdirp = os.path.join(featdir, "parts_input_cnsistent_clean")
        if( not os.path.exists(featdirp) ):
            os.makedirs(featdirp)

            # cnsistent acceps or the unique copy number value or the pair of minor & major copy number value, not both
            chunk = 20
            parts = {}
            path = os.path.join(featdir, "input_cnsistent.tsv")
            idx = 1
            ant = ""
            i = 0
            f = open(path, 'r')
            for line in f:
                l = line.replace('\n','').split('\t')
                line = '\t'.join(l[:4]+l[-2:])

                if(i==0):
                    header = line
                else:
                    if( (l[-1] not in ['nan', '2.0'] )and (l[-2] not in ['nan', '2.0']) ):
                        sample = line.split('\t')[0].split('_')[1]
                        if( sample != ant ):
                            if( (len(parts) % chunk == 0) and (len(parts) > 0) ) :
                                lines = list( map( lambda x: '\n'.join(x), parts.values() ))
                                lines = '\n'.join(lines)
                                opath = os.path.join(featdirp, "input_cnsistent_part%i.tsv" %(idx) )
                                g = open(opath, 'w')
                                g.write(header+"\n")
                                g.write(lines+"\n")
                                g.close()
                                
                                parts = {}
                                idx += 1
                            ant = sample


                        if(sample not in parts):
                            parts[sample] = []
                        parts[sample].append(line)
                i += 1

            if( len(parts) > 0):
                lines = list( map( lambda x: '\n'.join(x), parts.values() ))
                lines = '\n'.join(lines)
                path = os.path.join(featdir, "input_cnsistent_part%i.tsv" %(idx) )
                g = open(opath, 'w')
                g.write(header+"\n")
                g.write(lines+"\n")
                g.close()
            f.close()

    def _prepare_imputation_cnsistent(self, featdir):
        print("In imputation")
        indirp = os.path.join(featdir, "parts_input_cnsistent")
        indirp = os.path.join(featdir, "parts_input_cnsistent_clean")

        featdirp = os.path.join(featdir, "imputation_cnsistent")
        if( not os.path.exists(featdirp) ):
            os.makedirs(featdirp)
        
            print("In imputation")
            idx = 1
            for f in tqdm(os.listdir(indirp)):
                fin = os.path.join(indirp, f)
                fout = os.path.join(featdirp, "cns_imputed_part%i.tsv" %(idx) )
                if( not os.path.exists(fout) ):
                    print("\t", f)
                    cmd = "cns impute --method extend --assembly hg38 --out %s %s" %( fout, fin )
                    #os.system( cmd )
                    cmd = cmd.split(" ")
                    out = run( cmd, stdout=PIPE, stderr=PIPE )
                    if( out.returncode != 0 ):
                        raise Exception(out.stderr.decode())
                    else:
                        pass
                idx += 1

        # cns impute --method extend --assembly hg38 --out 'cns_test_out.tsv' /var/www/html/gdc_explorer_app/data_processed/TCGA-COAD/copy-number-variation/features/basic/parts_input_cnsistent/input_cnsistent_part1.tsv

    def _compute_cna_features_cnsistent(self, featdir):
        feats = ['coverage','ploidy','breakage']
        indirp = os.path.join(featdir, "imputation_cnsistent")
        
        #indirp = os.path.join(featdir, "parts_input_cnsistent_clean")
        featdirp = os.path.join(featdir, "results_cnsistent")
        if( not os.path.exists(featdirp) ):
            os.makedirs(featdirp)

            print("In extraction")
            idx = 1
            for f in tqdm(os.listdir(indirp)):
                fin = os.path.join(indirp, f)
                for fc in feats:
                    fout = os.path.join(featdirp, "cns_%s_part%i.tsv" %(fc, idx) )
                    if( not os.path.exists(fout) ):
                        cmd = "cns %s --assembly hg38 --out %s %s" %( fc, fout, fin )
                        #os.system( cmd )
                        cmd = cmd.split(" ")
                        out = run( cmd, stdout=PIPE, stderr=PIPE )
                        if( out.returncode != 0 ):
                            raise Exception(out.stderr.decode())
                        else:
                            pass
                idx += 1

    def _fuse_cnsistent_features(self, meta, featdir):
        # check if cnsistent sex inferred matches the patient metadata gender
        # check correlation between breakage and patient age

        indirp = os.path.join(featdir, "results_cnsistent")
        fused = os.path.join(featdir, "fused_features.tsv")
        if( not os.path.exists(fused) ):
            print("In fusion")
            dat = {}
            check_sex = {}
            check_breakage = { 'age': {}, 'breaks': {} }
            nsamples = 0
            for f in tqdm(os.listdir(indirp)):
                if( 'coverage' not in f ):
                    fin = os.path.join(indirp, f)
                    df = pd.read_csv(fin, sep='\t')
                    header = list(df.columns)[2:]
                    for i in df.index:
                        sample = df.loc[i, 'sample_id']
                        cond, case, uuid = sample.split('_')

                        if( not sample in dat ):
                            dat[sample] = { 'y': 0, 'y_vital': 0 }
                            if(cond == 'Tumor'):
                                dat[sample]['y'] = 1

                        if(case in meta):
                            check_sex[case] = 0
                            sex = df.loc[i, 'sex']
                            gender = meta[case]['gender'].lower()
                            age = int(meta[case]['age_at_initial_pathologic_diagnosis'])
                            check_breakage['age'][case] = age
                            if( (sex == 'xx' and gender == 'female') or (sex == 'xy' and gender == 'male') ):
                                check_sex[case] = 1

                            if( 'vital_status' in meta[case] ):
                                vital_status = meta[case]['vital_status']
                                if(vital_status == 'Dead'):
                                    dat[sample]['y_vital'] = 1
                        j = 0
                        for v in df.iloc[i, 2:]:
                            h = header[j]
                            dat[sample][h] = v
                            if( (case in meta) and (h == 'breaks_total_cn_all') ):
                                check_breakage['breaks'][case] = v
                            j += 1
            
            header = '\t'.join( ['sample_id'] + list(dat[sample]) )
            f = open(fused, 'w')
            f.write(header + '\n')
            for k in dat:
                values = [ str(x) for x in ([k] + list( dat[k].values() )) ]
                values = '\t'.join( values )
                f.write(values + '\n')
            f.close()

            # Results check sex
            nsamples = len(check_sex)
            print('cases: ', nsamples)
            matches = sum(check_sex.values())
            pc = (matches / nsamples)*100
            print('Match sex metadata x sex inferred: ', matches, ' - ', pc, '%')

            # Results correlation breakage
            breaks = list(check_breakage['breaks'].values())
            age = list(check_breakage['age'].values())
            print(len(breaks), len(age))
            corr, p_value = pearsonr(breaks, age)

            print('Correlation age x breaks: correlation - ', corr, ' | p-value - ', p_value)

    def extract_features_for_cnsistent(self, odir, fsodir, map_uuid_treat, map_case, meta):
        featdir = os.path.join(odir, "features", "cnsistent")
        if( not os.path.exists(featdir) ):
            os.makedirs(featdir)
        
        df = self._prepare_input_cnsistent( fsodir, map_uuid_treat, map_case, featdir)
        self._partitionate_input_cnsistent(featdir)
        #self._partitionate_input_cnsistent_clean_missing_neutral(featdir)
        self._prepare_imputation_cnsistent(featdir)
        self._compute_cna_features_cnsistent(featdir)
        self._fuse_cnsistent_features( meta, featdir)

    def test_feature_selection_rfecv(self, odir, path):
        # https://www.geeksforgeeks.org/machine-learning/recursive-feature-elimination-with-cross-validation-in-scikit-learn/
        
        df = pd.read_csv(path, sep='\t')
        x_cols = list(filter( lambda x: (not x.startswith('y') and not x.startswith('sample')), df.columns ))
        x = df.loc[:, x_cols].values
        y = df['y'].values

        classifiers = { 
            'decision_tree': DecisionTreeClassifier( random_state=0),
            'adaboost': AdaBoostClassifier(n_estimators=100, random_state=0)
        }
        print(df.groupby('y_vital').count())

        opath = os.path.join(odir, "results_rfecv.tsv")

        header = ['classifier', 'metric', 'value', 'selected_features', 'reduction_fraction']
        lines = [ header ]
        for clf in classifiers:
            for sel in [True, False]:
                print(clf, sel)
                estimator = classifiers[clf]
                
                selected = 'all'
                x_tmp = x
                rfrac = 0
                if(sel):
                    selector = RFECV(estimator, cv=10)
                    selector = selector.fit(x, y)
                    x_tmp = x[:, selector.support_]

                    # Print the optimal number of features
                    print("Optimal number of features: %d" % selector.n_features_)
                    rfrac = 1 - (selector.n_features_ / len(x_cols))

                    # Print the selected features
                    selected = ','.join(np.array(x_cols)[selector.support_])
                    print("Selected features: %s" %( selected ) )

                estimator.fit(x_tmp, y)
                
                metrics = ['f1','accuracy','precision','recall']
                for mt in metrics:
                    scores = cross_val_score(estimator, x, y, cv=10)
                    value = st.mean(scores)
                    lines.append( [clf, mt, value, selected, rfrac] )

        lines = list( map( lambda x: '\t'.join( [ str(y) for y in x ] ), lines ) )
        f = open( opath, 'w' )
        f.write( '\n'.join(lines) + '\n' )
        f.close()

    def _get_gene_civicdb_labels_for_diseases(self, odir, project):
        # the negative examples could be the random assignment of genes with low association score in opentarget

        stats = { 'genes': set(), 'evtype': set(), 'significance': set() }
        opath = os.path.join(odir, "y_civic.json")
        finaldb = {}
        ys = {}
        if( True or not os.path.exists(opath) ):
            from handle_geneset_enrichment import HandleEnrichment
            hge = HandleEnrichment()
            info = hge._load_civic_genes()
            mapp, revmapp = hge._get_map_project_diseaseCivic()
            diseases_filter = mapp[project]
            
            filtered_items = list( filter( lambda x: len( set(x[1]['disease']).intersection( set(diseases_filter) ) ) > 0, list(info.items()) ))
            for it in filtered_items:
                sep = ' AND '
                if( ' OR ' in it[0] ):
                    sep = ' OR '
                genes = it[0].split(sep)
                for gs in genes:
                    its = gs.split(' ')
                    ext = ''
                    gs = its[0]
                    if( len(its) > 1 ):
                        ext = its[1]
                    flag_fusion = ( '::' in gs )
                    for g in gs.split('::'): # fusion
                        if( not g in finaldb):
                            finaldb[g] = []
                            ys[g] = { 'evtype': set(), "significance": set() }
                        dat = it[1]
                        dat['comment'] = ext
                        dat['flag_fusion'] = flag_fusion
                        evtype = dat['evtype']
                        significance = dat['significance']
                        finaldb[g].append( dat )
                        ys[g]["evtype"].add(evtype)
                        ys[g]["significance"].add(significance)

                        stats['genes'].add(g)
                        stats['evtype'].add(evtype)
                        stats['significance'].add(significance)

            tmp = {}
            for g in ys:
                tmp[g] = { "evtype": list(ys[g]['evtype']), "significance": list(ys[g]['significance']) }
            dat = { "complete": finaldb, 'ys': tmp }
            json.dump( dat, open(opath, 'w') )
        else:
            dat = json.load( open(opath, 'r') )
            finaldb = dat['complete']
            ys = dat['ys']    
            for g in ys:
                stats['genes'].add(g)
                stats['evtype'].update( ys[g]['evtype'] )
                stats['significance'].update( ys[g]['significance'] )

        print('Stats')
        for k in stats:
            print(k, len(stats[k]), list(stats[k])[:10] )

        return finaldb, ys

    def _get_exp(self):
        project = 'TCGA-COAD'
        datcat = 'copy number variation'
        dattype = 'Gene Level Copy Number'

        # Gathering the directory for the project under CNA aspect
        odir, fsodir, file_list = self.proc.get_case_files_by_data_category(project, datcat)

        # Selecting only the files whose case has a tumor and a normal sample
        ok_cases, map_uuid_treat = self.proc._select_cases_with_normalAndTumor(odir, data_type = dattype)
        print( "Cases with normal x tumor", len(ok_cases) )
        # Retrieving raw files
        for uuid in map_uuid_treat:
            self.proc._get_file_by_uuid( fsodir, uuid)

        # Getting case annotation mappings
        map_case = self.proc._get_map_case_file(odir)
        meta = self.proc._get_cases_metadata(project)
        
        # Normalizing raw data with standard column names in the header
        #self._check_absolute_copy_number_value(fsodir, map_uuid_treat)

        # Trying out vanilla approach to transform the features
        #self.extract_common_features(odir, fsodir, map_case, map_uuid_treat, meta)
        
        # Test the possibility of prioritizing features according to the p-value comparing the values of a feature for all cases, in normal and tumor groups
        #self.test_significance_cnv_values(odir)

        # Trying the feature exraction approach based on cnsistent
        '''
        self.extract_features_for_cnsistent( odir, fsodir, map_uuid_treat, map_case, meta)
        featdir = os.path.join(odir, "features", "cnsistent")
        path = os.path.join(featdir, "fused_features.tsv")
        self.test_feature_selection_rfecv(featdir, path)
        '''

        self._get_gene_civicdb_labels_for_diseases(odir, project)


    def _filter_basic_data_type(self):
        # check coverage of cases per data_type
        project = 'TCGA-COAD'
        datcat = 'copy number variation'
        dattype = 'Gene Level Copy Number'
        odir, fsodir, file_list = self.proc.get_case_files_by_data_category(project, datcat)
        path = os.path.join(odir, "files_metadata.tsv")
        df = pd.read_csv(path, sep='\t')
        n_cases = len(df['cases.0.case_id'].unique())
        n_files = len(df['file_id'].unique())
        print('Total cases retrieved in files metadata:', n_cases)
        print('Total files retrieved in files metadata:', n_files)
        print(df.data_type.unique())

        dt = {}
        for i in df.index:
            case = df.loc[i, 'cases.0.case_id']
            dtype = df.loc[i, 'data_type']
            uuid = df.loc[i, 'file_id']
            condition = df.loc[i, 'cases.0.samples.0.tissue_type']
            if(not dtype in dt):
                dt[dtype] = { 'per_case': {}, 'per_condition': {} }

            if(not case in dt[dtype]['per_case'] ):
                dt[dtype]['per_case'][case] = {}
            if(not condition in dt[dtype]['per_case'][case]):
                dt[dtype]['per_case'][case][condition] = list()
            dt[dtype]['per_case'][case][condition].append(uuid)

            if(not condition in dt[dtype]['per_condition'] ):
                dt[dtype]['per_condition'][condition] = {}
            if(not case in dt[dtype]['per_condition'][condition]):
                dt[dtype]['per_condition'][condition][case] = list()
            dt[dtype]['per_condition'][condition][case].append(uuid)

        for dtype in dt:
            print('-------->', dtype)
            print("\tNumber of cases for ", dtype, ':', len(dt[dtype]['per_case']) )
            
            print("\tNumber of cases with files for Normal tissues for ", dtype, ':', len(dt[dtype]['per_condition']['Normal']) )
            print("\tNumber of cases with files for Tumor tissues for ", dtype, ':', len(dt[dtype]['per_condition']['Tumor']) )

            _cases = dt[dtype]['per_condition']['Normal'].keys()
            nfiles_normal = sum( [ len(dt[dtype]['per_condition']['Normal'][x]) for x in _cases ] )
            _cases = dt[dtype]['per_condition']['Tumor'].keys()
            nfiles_tumor = sum( [ len(dt[dtype]['per_condition']['Tumor'][x]) for x in _cases ] )
            print("\tNumber of files for Normal tissues for ", dtype, ':', nfiles_normal )
            print("\tNumber of files for Tumor tissues for ", dtype, ':', nfiles_tumor)

        print('Filtering uuids by gene level')
        ok_cases, map_uuid_treat = self.proc._select_cases_with_normalAndTumor(odir, data_type = dattype)
        print('\tNormal', len( list( filter( lambda x: x=='Normal', map_uuid_treat.values() )) ) )
        print('\tTumor', len( list( filter( lambda x: x=='Tumor', map_uuid_treat.values() )) ) )
        print('Total files by gene level', len(map_uuid_treat) )

        """
        Results:
            Total cases retrieved in files metadata: 461
            Total files retrieved in files metadata: 4655
            ['Gene Level Copy Number' 'Masked Copy Number Segment'
             'Copy Number Segment' 'Allele-specific Copy Number Segment']
            --------> Gene Level Copy Number
                Number of cases for  Gene Level Copy Number : 459
                Number of cases with files for Normal tissues for  Gene Level Copy Number : 215
                Number of cases with files for Tumor tissues for  Gene Level Copy Number : 453
                Number of files for Normal tissues for  Gene Level Copy Number : 444
                Number of files for Tumor tissues for  Gene Level Copy Number : 957
            --------> Masked Copy Number Segment
                Number of cases for  Masked Copy Number Segment : 460
                Number of cases with files for Normal tissues for  Masked Copy Number Segment : 416
                Number of cases with files for Tumor tissues for  Masked Copy Number Segment : 452
                Number of files for Normal tissues for  Masked Copy Number Segment : 470
                Number of files for Tumor tissues for  Masked Copy Number Segment : 506
            --------> Copy Number Segment
                Number of cases for  Copy Number Segment : 461
                Number of cases with files for Normal tissues for  Copy Number Segment : 427
                Number of cases with files for Tumor tissues for  Copy Number Segment : 455
                Number of files for Normal tissues for  Copy Number Segment : 638
                Number of files for Tumor tissues for  Copy Number Segment : 686
            --------> Allele-specific Copy Number Segment
                Number of cases for  Allele-specific Copy Number Segment : 454
                Number of cases with files for Normal tissues for  Allele-specific Copy Number Segment : 215
                Number of cases with files for Tumor tissues for  Allele-specific Copy Number Segment : 249
                Number of files for Normal tissues for  Allele-specific Copy Number Segment : 444
                Number of files for Tumor tissues for  Allele-specific Copy Number Segment : 510
            Filtering uuids by gene level
                Normal 433
                Tumor 229
        """

    def run(self):
        self._get_exp()
        
        featdir = '/var/www/html/gdc_explorer_app/data_processed/TCGA-COAD/copy-number-variation/features/basic/'
        #self._reconstruct_feat_cna_json(featdir)
        #self.analysis_cna_events(featdir)
        #self._filter_basic_data_type()

if( __name__ == "__main__" ):
    gtf_hs = '/mnt/yasdata/home/yasmmin/Dropbox/portfolio_2025/gdc_explorer_app_web/gdc_exploration_project/copy_number_variation/Homo_sapiens.GRCh38.115.gtf'
    out = '/mnt/yasdata/home/yasmmin/Dropbox/portfolio_2025/gdc_explorer_app_web/gdc_exploration_project/copy_number_variation/out'
    
    o = HandleCNAAnalysis( out, gtf = gtf_hs )
    o.run()