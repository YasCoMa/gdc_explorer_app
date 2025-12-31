import os
import sys
import json
import pickle
import requests

import numpy as np
import pandas as pd
import gseapy as gp

from tqdm import tqdm
from Bio import Entrez
import xml.etree.ElementTree as ET

from pydeseq2.dds import DeseqDataSet
from pydeseq2.default_inference import DefaultInference
from pydeseq2.ds import DeseqStats

class HandleEnrichment:
    # Library reference: https://gseapy.readthedocs.io/en/latest/gseapy_example.html#GSEA-Example

    def __init__(self, genes_prey=None, fout = '../data_processed'):
        Entrez.email = 'ycfrenchgirl2@gmail.com'
        Entrez.api_key="4543094c8a41e6aecf9a1431bff42cfac209"

        self.out = fout
        if( not os.path.exists(self.out) ):
            os.makedirs( fout )

        libs = { "h": "hallmark", "c1": "positional", "c2": "curated", "c3": "regulatory", "c4": "computational", "c5": "ontology", "c6": "oncogenic", "c7": "immunologic", "c8": "cell types" }
        self.prey = genes_prey


    def _get_significance_dbsnp(self, rsid):
        
        rsid = rsid.replace('rs', '')
        fetch = Entrez.efetch(db='snp',resetmode='xml',id=rsid,rettype='full')
        t = fetch.read()
        tree = ET.fromstring(t)
        nodes = list(tree.iter())
        try:
            significance = list( filter(lambda x: x.tag.lower().endswith('significance'), nodes) )[0].text
        except:
            significance = None

        return significance

    def _get_info_clinvar(self, rsid):
        handle = Entrez.esearch(db='clinvar', sort='relevance', term=rsid, retmax=100)
        res = Entrez.read(handle)
        cvids = res['IdList']

        dat = { 'snp_id': rsid, 'clinvar_assoc': {} }
        for c in cvids:
            diseases = set()
            url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=clinvar&rettype=vcv&is_variationid&id=%s" %(c)
            r = requests.get( url )
            tree = ET.fromstring(r.text)
            nodes = list(tree.iter())

            cls_sections = list( filter(lambda x: x.tag.lower() == 'classifications', nodes) )
            print( c, len(cls_sections) )
            aux = []
            for sec in cls_sections:

                inner_nodes = list(sec.iter())
                significance = list( filter(lambda x: x.tag.lower() == 'description', inner_nodes) )[0].text
                obj = {'significance': significance, 'diseases': set() }

                traits = list( filter(lambda x: x.tag.lower() =='trait' and x.attrib['Type']=='Disease', inner_nodes) )
                
                for tr in traits:
                    tr_nodes = list(tr.iter())
                    names = list( filter(lambda x: x.tag.lower() == 'elementvalue' and x.attrib['Type']=='Preferred', tr_nodes) )
                    names = list( map( lambda x: x.text, names ) )
                    diseases.update(names)
                    obj['diseases'].update(names)
                aux.append(obj)

            dat['clinvar_assoc'][c] = { "assoc_significance": aux, "all_diseases": diseases }
        return dat

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

    
    def test_differential_expression(self, outdir, ide, metadata, counts_df, cutoff_log2fc=2, cutoff_padj=0.001):
        '''
        outdir = './'
        indir = '/var/www/html/gdcexplorer_app/data_processed/TCGA-BLCA_transcriptome-profiling/'

        # Loading
        cpath = os.path.join(indir, 'deseq_table_counts.tsv')
        counts_df = pd.read_csv( cpath, sep='\t', index_col=0)

        mpath = os.path.join(indir, 'deseq_table_meta.tsv')
        metadata = pd.read_csv( mpath, sep='\t', index_col=0)
        '''

        # Preprocessing data
        inference = DefaultInference(n_cpus=6)
        dds = DeseqDataSet( counts=counts_df, metadata=metadata, design="~condition", refit_cooks=True, inference=inference, )

        dds.fit_size_factors()
        # dds.obs["size_factors"]

        dds.fit_genewise_dispersions()
        # dds.var["genewise_dispersions"]

        dds.fit_dispersion_trend()
        # dds.uns["trend_coeffs"]
        # dds.var["fitted_dispersions"]

        dds.fit_dispersion_prior()
        # print( f"logres_prior={dds.uns['_squared_logres']}, sigma_prior={dds.uns['prior_disp_var']}" )

        dds.fit_MAP_dispersions()
        # dds.var["MAP_dispersions"]
        # dds.var["dispersions"]

        dds.fit_LFC()
        # dds.varm["LFC"]

        dds.calculate_cooks()
        dds.refit()

        # Performing stats significance tests
        ds = DeseqStats(dds, contrast=np.array([0, 1]), alpha=0.05, cooks_filter=True, independent_filter=True,)
        ds.run_wald_test()
        # ds.p_values

        if( ds.cooks_filter ):
            ds._cooks_filtering()
        #ds.p_values

        if ds.independent_filter:
            ds._independent_filtering()
        else:
            ds._p_value_adjustment()
        # ds.padj

        ds.summary()
        de = ds.results_df

        opath = os.path.join(outdir, '%s_deg_result.pkl' %(ide) )
        pickle.dump( ds, open(opath, 'wb') )

        # Post processing
        downs = de[ (de.log2FoldChange <= (-1*cutoff_log2fc) ) & ( de.padj <= 0.001 ) ]
        aux = {}
        keys = list(downs.index)
        for k in keys:
            aux[k] = { 'log2fc': downs.loc[k, 'log2FoldChange'], 'padj': downs.loc[k, 'padj'] }
        downs = aux
        
        ups = de[ (de.log2FoldChange >= cutoff_log2fc ) & ( de.padj <= 0.001 ) ]
        aux = {}
        keys = list(ups.index)
        for k in keys:
            aux[k] = { 'log2fc': ups.loc[k, 'log2FoldChange'], 'padj': ups.loc[k, 'padj'] }
        ups = aux

        result = { 'up': ups, 'downs': downs }
        oname = "%s_result_log2fc-%s_padj-%s.json" %( ide, str(cutoff_log2fc), str(cutoff_padj).split('.')[1] )
        opath = os.path.join( outdir, oname )
        json.dump( result, open(opath, 'w') )

    def _test_proportion_demovar(self, mdf, dim, cutoff=0.7):
        flag = True

        n = len(mdf)
        max_part = 0
        subgroups = mdf[dim].unique()
        for s in subgroups:
            aux = mdf[ mdf[dim] == s ]
            sp = len(aux)
            portion = sp/n
            if(portion > max_part):
                max_part = portion
        
        if(max_part > cutoff):
            flag = False

        return flag        

    def perform_degs_analysis_simulation(self, project):
        print('------->', project)

        datcat = "transcriptome profiling"
        basename = "%s_%s" %(project, datcat.replace(" ", "-"))
        indir = os.path.join(self.out, "%s" %(basename) )
        outdir = os.path.join(self.out, "%s" %(basename), 'deg_analysis' )
        if( not os.path.isdir( outdir ) ):
            os.makedirs( outdir )

        cpath = os.path.join(indir, 'deseq_table_counts.tsv')
        if( os.path.exists(cpath) ):
            counts_df = pd.read_csv( cpath, sep='\t', index_col=0)

            mpath = os.path.join(indir, 'deseq_table_meta.tsv')
            metadata = pd.read_csv( mpath, sep='\t', index_col=0)

            ide = "by_all"
            self.test_differential_expression(outdir, ide, metadata, counts_df)
            print( '\t', 'all', len(metadata) )

            cols_stratification = ['race','gender', 'ethnicity']
            for c in cols_stratification:
                flag = self._test_proportion_demovar(metadata, c)
                if(flag):
                    k = "by_%s" %(c)
                    aux_outdir = os.path.join( outdir, k )
                    if( not os.path.isdir( aux_outdir ) ):
                        os.makedirs( aux_outdir )

                    subgroups = metadata[c].unique()
                    for s in subgroups:
                        ide = "by_%s-group_%s_" %(c, s)
                        meta_aux = metadata[ metadata[c] == s ]
                        print( '\t', c, s, len(meta_aux) )

                        samples = list( meta_aux.index )
                        counts_aux = counts_df.iloc[samples, :]
                        self.test_differential_expression( aux_outdir, ide, meta_aux, counts_aux)
            print('\n')

    def run(self):
        p = 'TCGA-ACC'
        #self.enrich_exclusive_mutated_genes(p)

        #d = self._get_info_clinvar('rs334')
        #print(d)

        projects = [ "TCGA-ACC",  "TCGA-BLCA",  "TCGA-BRCA",  "TCGA-CESC",  "TCGA-CHOL",  "TCGA-COAD",  "TCGA-DLBC",  "TCGA-ESCA",  "TCGA-GBM",  "TCGA-HNSC",  "TCGA-KICH",  "TCGA-KIRC",  "TCGA-KIRP",  "TCGA-LAML",  "TCGA-LGG",  "TCGA-LIHC",  "TCGA-LUAD",  "TCGA-LUSC",  "TCGA-MESO",  "TCGA-OV",  "TCGA-PAAD",  "TCGA-PCPG",  "TCGA-PRAD",  "TCGA-READ",  "TCGA-SARC",  "TCGA-SKCM",  "TCGA-STAD",  "TCGA-TGCT",  "TCGA-THCA",  "TCGA-THYM",  "TCGA-UCEC",  "TCGA-UCS",  "TCGA-UVM" ]
        for p in tqdm(projects):
            self.perform_degs_analysis_simulation(p)

if( __name__ == "__main__" ):
    o = HandleEnrichment()
    o.run()
