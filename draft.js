
mapp = await obj_cov.get_methylation_mapping()

p = 'TCGA-KIRC'
datcat = 'dna methylation'
dat = await obj_cov.get_case_files_by_data_category(p, datcat)
uuids = dat.filter( e => e.file_size <= obj_cov.file_size_limit ).map( e => e.file_id )

// mapping cpg ids to genes and chromosomes: HumanMethylation450 v1.2 Manifest File (CSV Format) in https://support.illumina.com/downloads/infinium_humanmethylation450_product_files.html

// classifying files per tissue and tumor
	y = {  }

	// stratifying the files according to the subgroups
	filesdf = { "race": {}, "gender": {}, "ethnicity": {} }

	// How many had recurrence?
	recurrents = []

	// In the cases that there was a tumor, what was the mean diagnosis age per subgroup?
	samples = {}

	dat.forEach( el => {
		let uuid = el.file_id;
		let size = el.file_size;
		let platform = el.platform;

		try{
			if( size <= obj_cov.file_size_limit ) {
				let sample_id = el.cases[0].submitter_id;
				let oaux = { 'uuid': uuid, 'submitter_id': sample_id, "size": size, "platform": platform };
				let sobj = { "sample_id": sample_id };

				let tissue = el.cases[0].samples[0].tissue_type;
				let stage = el.cases[0].samples[0].tumor_descriptor;
				y[uuid] = { "tissue": tissue, "stage": stage };

				let races = Object.keys( filesdf['race'] );
				let dr = el.cases[0].demographic.race;
				sobj["race"] = dr;
				if( ! races.includes(dr) ){
					filesdf['race'][dr] = [];
				}
				filesdf['race'][dr].push( oaux );

				let genders = Object.keys( filesdf['gender'] );
				dr = el.cases[0].demographic.gender;
				sobj["gender"] = dr;
				if( ! genders.includes(dr) ){
					filesdf['gender'][dr] = [];
				}
				filesdf['gender'][dr].push( oaux );

				let ethncities = Object.keys( filesdf['ethnicity'] );
				dr = el.cases[0].demographic.ethnicity;
				sobj["ethnicity"] = dr;
				if( ! ethncities.includes(dr) ){
					filesdf['ethnicity'][dr] = [];
				}
				filesdf['ethnicity'][dr].push( oaux );

				samples[uuid] = sobj;

				let diags = el.cases[0].diagnoses;
				for( let dg of diags ){
					console.log(dg)
					if( dg.days_to_recurrence ){
						dg.age_at_diagnosis = parseInt( dg.age_at_diagnosis / 365 );
						dg.days_to_recurrence = parseInt( dg.days_to_recurrence / 365 );
						dg.days_to_last_follow_up = parseInt( dg.days_to_last_follow_up / 365 );
						let obj = dg;
						obj['uuid'] = uuid;
						obj['submitter_id'] = sample_id;
						recurrents.push( obj );
					}
				}
			}
		}
		catch(e){
			console.log(el)
			console.log('error', e)
		}

	} )

// Filter samples whose file_size is under 20Mb

// Checking tumor / normal samples distribution
	tumor = Object.keys(y).filter( e => y[e].tissue == 'Tumor' ).length
	normal =  Object.keys(y).filter( e => y[e].tissue == 'Normal' ).length


dfs = await obj_cov.retrieve_process_methylation_files( uuids )

f = await obj_cov._get_file_by_uuid("daa01486-58c6-4fc3-aab5-a0c4680ee10f") // example of file with 765 Kb