import os

src = "data_processed"
dest = "data_release"

def transfer():
	targets = [
		"all_map.json",
		"_project_/clinical/data_cases.json",
		"_project_/clinical/survival_probs.json",
		"_project_/cases/cases_metadata.json",
		"_project_/simple-nucleotide-variation/files_metadata.tsv",
		"_project_/simple-nucleotide-variation/general_snprs_info.json",
		"_project_/simple-nucleotide-variation/by_race_table_summary.tsv",
		"_project_/simple-nucleotide-variation/by_gender_table_summary.tsv",
		"_project_/simple-nucleotide-variation/by_ethnicity_table_summary.tsv",
		"_project_/simple-nucleotide-variation/pgx_data.json",
		"_project_/simple-nucleotide-variation/clingx_filtered_annotation_info_for_cases.json"
	]
	targets = [
		"_project_/simple-nucleotide-variation/clingx_filtered_annotation_info_for_cases.json"
	]

	projects = list( filter( lambda x: x.startswith("TCGA"), os.listdir(src) ) )
	for p in projects:
		for t in targets:
			t = t.replace("_project_", p)
			spath = os.path.join(src, t)

			ddir = '/'.join( t.split('/')[:-1] )
			dfile = t.split('/')[-1]
			ddir = os.path.join( dest, ddir )
			if( not os.path.exists(ddir) ):
				os.makedirs(ddir)
			dpath = os.path.join(ddir, dfile)

			os.system( f"cp {spath} {dpath}" )

transfer()