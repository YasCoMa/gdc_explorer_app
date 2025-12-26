import os
import sys

import numpy as np
from gseapy import Msigdb

class HandleEnrichment:
	# Library reference: https://gseapy.readthedocs.io/en/latest/gseapy_example.html#GSEA-Example

	def __init__(self, genes_prey):
		self.prey = genes_prey


	def init_msigdb(self):
		# Retrieves last release of the hallmark gene sets from msigdb
		msig = Msigdb()
		gmt = msig.get_gmt(category='h.all', dbver="2025.1.Hs")
		'''
		Some geneset important for cancer omics analysis are not in enrichr or in this msig collection, but they can be downloaded from gsea site:
		https://www.gsea-msigdb.org/gsea/msigdb/cards/KEGG_HOMOLOGOUS_RECOMBINATION
		https://www.gsea-msigdb.org/gsea/msigdb/cards/GOBP_HOMOLOGOUS_RECOMBINATION

		then you can pass as a path: gene_sets=["./tests/data/khr.gmt", ...]

		names = gp.get_library_name()
		'''

	
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
		ss = gp.ssgsea(data = gs, gene_sets=gmt, outdir=None, sample_norm_method='rank', no_plot=True)
        df = ss.res2d

        return df

    def check_genes_covariance(self, inn)
    	# covariance reference: https://www.sciencedirect.com/science/article/pii/S0092867418307232#fig3
    	#	they plotted the distribution of variance across the groups being compared (normal, tumor)
    	#	this article was the first to describe and define phenotypic volume	

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
