import os
import re
import json
import requests

import pandas as pd

"""
Goal: Provide a tool to treat and augment tcga/gdc data from the bioStudies platform

todo: provide statistics of class cna
"""

class IntegrationCNABioStudies:
    def __init__(self, fout, gtf = None):
        self.out = fout
        if( not os.path.exists(self.out) ):
            os.makedirs( self.out )

        self.gtf_map = None
        if( (gtf is not None) and (os.path.exists(gtf)) ):
            self.gtf_map = self.process_gtf(gtf)

        self.processed = os.path.join( self.out, "processed" )
        if( not os.path.exists(self.processed) ):
            os.makedirs( self.processed )

        self.raw_files = os.path.join( self.out, "raw_files" )
        if( not os.path.exists(self.raw_files) ):
            os.makedirs( self.raw_files )
    
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

    def query_studies(self, cancer_type):
        studies = []

        page = 1
        pgsize = 100
        u = f"https://www.ebi.ac.uk/biostudies/api/v1/CancerModelsOrg/search?pageSize={pgsize}&sortBy=relevance&page={page}&title={cancer_type}&facet.cancermodel\
sorg.dataset_available=copy%20number%20alteration"
        rs = requests.get(u)
        d = rs.json()
        
        next_pages = []
        total = d['totalHits']
        if( total > 0 ):
            studies += list( map( lambda x: x['accession'], d['hits'] ))

            if( total > pgsize ):
                npages = total//pgsize
                if( total % pgsize > 0 ):
                    npages += 1
                next_pages = list( range(2, npages+1) ) 
                flag = True

        
        for page in next_pages:
            u = f"https://www.ebi.ac.uk/biostudies/api/v1/CancerModelsOrg/search?pageSize={pgsize}&sortBy=relevance&page={page}&title={cancer_type}&facet.cancermodel\
sorg.dataset_available=copy%20number%20alteration"
            rs = requests.get(u)
            d = rs.json()
            studies += list( map( lambda x: x['accession'], d['hits'] ))

        obj = { "cancer_type": cancer_type, "result_studies": studies }
        opath = os.path.join(self.out, "result_query.json")
        json.dump(obj, open(opath, 'w') )

    def __download_raw_file(self, accid, uuid, filepath):
        opath = os.path.join( self.raw_files, "%s_raw.tsv" %(uuid) )
        if( True or not os.path.exists(opath) ):
            parts = re.findall( r'(\d+)', accid)
            nid = parts[0]
            p1 = accid.replace( nid, '' )
            p2 = nid[-3:]
            while( len(p2) != 3 ):
                p2 = '0'+p2
            ext = filepath.split('/')[-1]
            u = "https://ftp.ebi.ac.uk/pub/databases/biostudies/%s/%s/%s/Files/%s/%s" %( p1, p2, accid, filepath, ext )
            rs = requests.get(u)
            txt = rs.text
            
            f = open(opath, 'w')
            f.write(txt)
            f.close()

    def __check_cna_datatype(self, attrs):
        flag = False
        obj = {  }
        for a in attrs:
            if(a['name'] == 'Sample ID'):
                obj['sample_id'] = a['value']
            if(a['name'] == 'Sample Type'):
                obj['sample_type'] = a['value']

            if(a['name'] == 'Data Type'):
                if(a['value'] == 'copy number alteration'):
                    flag = True

        return flag, obj

    def retrieve_study_information(self, accid):
        opath = os.path.join(self.out, "result_parsing_metadata_%s.json" %(accid) )
        if(not os.path.exists(opath) ):
            meta = { "patient_tumor_info": {}, "files": {} }
            u = "https://www.ebi.ac.uk/biostudies/api/v1/studies/%s" %(accid)
            rs = requests.get(u)
            d = rs.json()

            for s in d['section']['subsections']:
                if( 'type' in s ):
                    if( s['type'] == 'Patient / Tumour Metadata' ):
                        for attr in s['attributes']:
                            key = attr['name'].lower().replace(' ', '_')
                            value = attr['value'].lower().replace(' ', '_')
                            meta["patient_tumor_info"][key] = value

                if( 'files' in s ):
                    for f in s['files']:
                        attrs = f['attributes']
                        flag, obj = self.__check_cna_datatype( attrs)
                        if( flag ):
                            obj['path'] = f['path']
                            meta["files"][ obj['sample_id'] ] = obj
                            uuid = "%s@%s" %(accid, obj['sample_id'] )
                            self.__download_raw_file( accid, uuid, obj['path'])

            json.dump(meta, open(opath, 'w') )

    def __get_class_cnv(self, cnv):
        class_cnv = 'neutral'
        if( cnv >= 2.64 and cnv <= 3.36 ):
            class_cnv = 'gain'
        elif( cnv > 3.36 ):
            class_cnv = 'amp'
        elif( cnv >= 0.87 and cnv <= 1.32 ):
            class_cnv = 'loss'
        elif( cnv < 0.87 ):
            class_cnv = 'del'

        return class_cnv

    def get_absolute_copy_number_harmonize(self, uuid):
        # https://cnvkit.readthedocs.io/en/stable/pipeline.html#call
        opath = os.path.join( self.processed, "%s_processed.tsv" %(uuid) )
        if(not os.path.exists(opath) ):
            path = os.path.join( self.raw_files, "%s_raw.tsv" %(uuid) )
            print(path)
            toparse = {}
            dat = {}
            df = pd.read_csv(path, sep='\t')
            for i in df.index:
                gene = df.loc[i, 'hgnc_symbol']
                if( gene in self.gtf_map ):
                    cnv = df.loc[i, 'gistic_value']
                    log2r = df.loc[i, 'log2r_cna']
                    
                    gid = self.gtf_map[gene]['gid']
                    chrom = self.gtf_map[gene]['chrom']
                    start, end = self.gtf_map[gene]['coord']
                    dat[gene] = [ gid, gene, chrom, start, end, cnv, cnv, cnv ]
                    if( (str(cnv) == 'nan') and ( str(log2r) != 'nan' ) ):
                        toparse[gene] = [ gid, chrom, start, end, gene, log2r ]
                    

            info = {}
            if( len(toparse) > 0 ):
                header = '\t'.join(['chromosome','start','end', 'gene', 'log2'] )
                op = 'test_cnv.cns'
                f = open(op,'w')
                f.write(header + '\n')
                for g in toparse:
                    line = '\t'.join( [ str(x) for x in toparse[g][1:] ] ) + '\n'
                    f.write(line)
                f.close()

                os.system('cnvkit.py call %s -o tmp.cns' %(op) )
                ip = 'tmp.cns'
                df = pd.read_csv(ip, sep='\t')
                info = dict( zip( df.gene.values.tolist(), df.cn.values.tolist() ))

            header = '\t'.join(['gene_id', 'gene_name', 'chromosome', 'start', 'end', 'copy_number', 'min_copy_number', 'max_copy_number', 'class_cnv'])
            f = open( opath, 'w' )
            f.write(header + '\n')
            for g in dat:
                v, minv, maxv = dat[g][-3:]

                if(g in info):
                    v = info[g]
                    minv = info[g]
                    maxv = info[g]
                
                try:
                    v = int(v)

                    class_cnv = self.__get_class_cnv(v)
                    values = dat[g][:-3] + [v, minv, maxv, class_cnv]
                    values = '\t'.join( [ str(x) for x in values ] )
                    f.write(values + '\n')
                except:
                    pass
            f.close()

    def process_query_results(self):
        path = os.path.join(self.out, "result_query.json")
        d = json.load( open(path, 'r') )
        for accid in d["result_studies"]:
            self.retrieve_study_information(accid)
            
            path = os.path.join(self.out, "result_parsing_metadata_%s.json" %(accid) )
            meta = json.load( open(path, 'r') )
            for fsample in meta["files"]:
                uuid = "%s@%s" %(accid, fsample)
                self.get_absolute_copy_number_harmonize(uuid)

    def run(self):
        ctype = "colon"
        #self.query_studies(ctype)
        self.process_query_results()
        
        '''
        accid = 'S-CMO92'
        uuid = 'S-CMO92@SIDS00003'
        filepath = 'CMP/SIDM00823/molecular_data/SIDS00003_copy-number-alteration_Illumina-HiSeq-2000.tsv'
        self.__download_raw_file(accid, uuid, filepath)
        '''

if( __name__ == "__main__" ):
    gtf_hs = '/mnt/yasdata/home/yasmmin/Dropbox/portfolio_2025/gdc_explorer_app_web/gdc_exploration_project/copy_number_variation/Homo_sapiens.GRCh38.115.gtf'
    
    out = '/mnt/yasdata/home/yasmmin/Dropbox/portfolio_2025/gdc_explorer_app_web/gdc_exploration_project/copy_number_variation/biostudies'
    
    o = IntegrationCNABioStudies( out, gtf = gtf_hs )
    o.run()