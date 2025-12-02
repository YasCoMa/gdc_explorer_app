import os
import sys
import json

class DataWrangler:
	def __init__(self, fout = '../data_processed'):
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

	def run(self):
		option = sys.argv[1]
		print(option)
		#try:
		eval('self.run_'+option)()
		#except:
		#	print("Function not available")

if( __name__ == "__main__" ):
	o = DataWrangler()
	o.run()
