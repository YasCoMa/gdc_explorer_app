import os
import sys
import json

import numpy as np
import gseapy as gp

class HandleEnrichment:
    # Library reference: https://gseapy.readthedocs.io/en/latest/gseapy_example.html#GSEA-Example

    def __init__(self, genes_prey=None, fout = '../data_processed'):
        self.out = fout
        if( not os.path.exists(self.out) ):
            os.makedirs( fout )

        libs = { "h": "hallmark", "c1": "positional", "c2": "curated", "c3": "regulatory", "c4": "computational", "c5": "ontology", "c6": "oncogenic", "c7": "immunologic", "c8": "cell types" }
        self.prey = genes_prey


    def init_msigdb(self, subset):
        # Retrieves last release of the hallmark gene sets from msigdb
        msig = gp.Msigdb()
        gmt = msig.get_gmt(category='%s.all' %(subset), dbver="2025.1.Hs")

        '''
        types of genesets (https://www.gsea-msigdb.org/gsea/msigdb/human/genesets.jsp):
        h - hallmark - gmt = msig.get_gmt(category='h.all', dbver="2025.1.Hs")
        c1 - positional (by chromosome) - gmt = msig.get_gmt(category='c1.all', dbver="2025.1.Hs")
        c2 - curated - gmt = msig.get_gmt(category='c2.all', dbver="2025.1.Hs")
        c3 - regulatory - gmt = msig.get_gmt(category='c3.all', dbver="2025.1.Hs")
        c4 - computational, cancer modules - gmt = msig.get_gmt(category='c4.all', dbver="2025.1.Hs")
        c5 - ontology - gmt = msig.get_gmt(category='c5.all', dbver="2025.1.Hs")
        c6 - oncogenic - gmt = msig.get_gmt(category='c6.all', dbver="2025.1.Hs")
        c7 - immunologic - gmt = msig.get_gmt(category='c7.all', dbver="2025.1.Hs")
        c8 - cell types - gmt = msig.get_gmt(category='c8.all', dbver="2025.1.Hs")

        Some geneset important for cancer omics analysis are not in enrichr or in this msig collection, but they can be downloaded from gsea site:
        https://www.gsea-msigdb.org/gsea/msigdb/cards/KEGG_HOMOLOGOUS_RECOMBINATION
        https://www.gsea-msigdb.org/gsea/msigdb/cards/GOBP_HOMOLOGOUS_RECOMBINATION
        https://www.gsea-msigdb.org/gsea/msigdb/cards/GOBP_REGULATION_OF_HORMONE_LEVELS

        then you can pass as a path: gene_sets=["./tests/data/khr.gmt", ...]

        names = gp.get_library_name()
        list( filter( lambda x: x.lower().find("hormone") != -1, names ))
        '''
        return gmt
    
    def enrich_genes_pvalue(self, genes):
        enr = gp.enrichr(gene_list=genes, # or "./tests/data/gene_list.txt",
                 gene_sets=['MSigDB_Hallmark_2020','KEGG_2021_Human'],
                 organism='human', # don't forget to set organism to the one you desired! e.g. Yeast
                 outdir=None, # don't write to disk
                )
        df = enr.res2d

        return df

    def enrich_single_sample_compare_phenotype(self, gs, gmt_library):
        # gs = { 'female_significants': [], 'male_significants': ['gene1', 'gene2'] }
        # or a table in which the columns Gene and NAME are the genes significant values, and the rows are the samples
        ss = gp.ssgsea(data = gs, gene_sets=gmt_library, outdir=None, sample_norm_method='rank', no_plot=True)
        df = ss.res2d

        return df

    def check_genes_covariance(self, inn):
        # covariance reference: https://www.sciencedirect.com/science/article/pii/S0092867418307232#fig3
        #    they plotted the distribution of variance across the groups being compared (normal, tumor)
        #    this article was the first to describe and define phenotypic volume    

        # inn is a dictionary whose key is the gene identifier, and the values in the array value are the counts for the aspect being analyzed (types of variant). Ex.: { 'missense': { 'gene1': [1,2,3,4,5,6], 'gene2': [34,2,7,4,3,10] } } , the size of the vector may be the number of samples in a subgroup (female). The arrays will always have the same value because in case a gene does not have a specific type

        # How much the cov matrix changes among the subgroups and taking all population?
        dat = {}
        for k in inn:
            sdf = {}
            covmat = [ inn[k][g] for g in inn[k] ]
            cov = np.cov(cov_mat)
            for i, g1 in enumerate(inn[k]):
                sdf[g1] = {}
                for j, g2 in enumerate(inn[k]):
                    sdf[g1][g2] = cov[i,j]

            dat[k] = sdf


        return dat

    def _get_exclusive_genes(self, df, cols_stratification):
        feature = 'hugo_symbol'
        all_exclusive = set()
        exclusive_dat = {}
        
        # Getting exclusive genes
        for c in cols_stratification:
            k = "by_%s" %(c)
            subgroups = list( df[k] )
            aux = {}
            for sroot in subgroups:
                subtotal = set()
                for s in subgroups:
                    if( s != sroot):
                        subtotal.update( list( df[k][s][feature] ) )
                aux[sroot] = set( df[k][sroot][feature] ) - subtotal
                all_exclusive = all_exclusive.union(aux[sroot])
            exclusive_dat[k] = aux

        return exclusive_dat, all_exclusive

    def _build_gtc_matrix(self, cols_stratification, odir, exclusive_genes, all_genes, df):
        feature = 'hugo_symbol'
        
        tables = {}
        for c in cols_stratification:
            k = "by_%s" %(c)

            subgroups = list( df[k] )
            header = ["Gene", "NAME"] + subgroups
            header = '\t'.join(header)
            lines = [ header ]
            aux = {}
            for g in all_genes:
                aux[g] = {}
            for s in subgroups:
                arr = []
                for g in all_genes:
                    if g in exclusive_genes[k][s]:
                        arr = df[k][s][feature][g].values()
                        aux[g][s] = sum(arr)/len(arr) # test min or max
                        aux[g][s] = 1
                    else:
                        aux[g][s] = 0

            for gene in aux:
                values = list( aux[gene].values() )
                l = [gene, gene ] + values
                l = [ str(x) for x in l ]
                lines.append( '\t'.join( l ) )

            opath = os.path.join(odir, "%s-enrich_table_exclusive_genes.tsv" %(k) )
            f = open(opath, 'w')
            f.write( '\n'.join(lines)+'\n' )
            f.close()
            tables[k] = opath

        return tables

    def enrich_exclusive_mutated_genes(self, project):
        datcat = "simple nucleotide variation"
        basename = "%s_%s" %(project, datcat.replace(" ", "-"))
        odir = os.path.join(self.out, "%s" %(basename) )
        path = os.path.join(odir, "data_cases.json")
        df = json.load( open(path, 'r') )

        cols_stratification = ['race','gender', 'ethnicity']
        exclusive_dat, all_exclusive = self._get_exclusive_genes(df, cols_stratification)
        
        # Enriching
        target_libs = ["c2", "c5", "c7"]
        tables = self._build_gtc_matrix(cols_stratification, odir, exclusive_dat, all_exclusive, df)
        for l in target_libs:
        	gmt = self.init_msigdb(l)
        	for k in tables:
        		inpath = tables[k]
	        	df = self.enrich_single_sample_compare_phenotype(inpath, gmt)
	        	opath = os.path.join(odir, "%s_result_%s-enrich_table.tsv" %(l, k) )
	        	df.to_csv(opath, sep='\t', index=None)

    def run(self):
        p = 'TCGA-ACC'
        self.enrich_exclusive_mutated_genes(p)
        
        projects = [ "TCGA-ACC",  "TCGA-BLCA",  "TCGA-BRCA",  "TCGA-CESC",  "TCGA-CHOL",  "TCGA-COAD",  "TCGA-DLBC",  "TCGA-ESCA",  "TCGA-GBM",  "TCGA-HNSC",  "TCGA-KICH",  "TCGA-KIRC",  "TCGA-KIRP",  "TCGA-LAML",  "TCGA-LGG",  "TCGA-LIHC",  "TCGA-LUAD",  "TCGA-LUSC",  "TCGA-MESO",  "TCGA-OV",  "TCGA-PAAD",  "TCGA-PCPG",  "TCGA-PRAD",  "TCGA-READ",  "TCGA-SARC",  "TCGA-SKCM",  "TCGA-STAD",  "TCGA-TGCT",  "TCGA-THCA",  "TCGA-THYM",  "TCGA-UCEC",  "TCGA-UCS",  "TCGA-UVM" ]
        #for p in tqdm(projects):
            #self.parse_clinical_data(p)
            #self.test_survival_km(p, dc)
            #self._compress_files(p, dc, remove=True)
            #self.parse_mutationSnv_data(p)

if( __name__ == "__main__" ):
    o = HandleEnrichment()
    o.run()
