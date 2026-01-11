class GdcExplorerLib {

    constructor(){
        this.file_size_limit = 1000000;
        this.meta = { 'program': [], 'disease': [], 'primary_site': [] };
        this.pids = [];
        this.projects = [];
    	this.projects_summary = {};
    	this.processed_counts = {};
        this.methylation_cgids_ann = {};
    	this.formats_datCategory = { "biospecimen": ["svs", "jpeg 2000"], "clinical": ["bcr xml"], "copy number variation": ["tsv", "txt"], "dna methylation": ["txt"], "proteome profiling": ["tsv"], "simple nucleotide variation": ["maf"], "transcriptome profiling": ["tsv"] };
    	
    }
    // Coverage is given by the case_couunt in each entry of exp strategy or data category divided by the general case_count of the project
    /*
    The combined nuleotide variation, somatic structural variation and structural variation data types contain only controlled access files. The sequencing reads type not only it is controlled but it needs a whole pipeline to be ready to use in downstream analysis.
    
    */
    
    get_metadata_options(){
        let programs = [];
        let diseases = [];
        let sites = [];
        for( let p of this.projects ){
            if( p.project_id.includes('-') ){
                let prg = p.project_id.split('-')[0]
                if( ! programs.includes(prg) ){
                    programs.push(prg);
                }
            }

            for( let d of p.disease_type ){
                if( ! diseases.includes(d) ){
                    diseases.push(d);
                }
            }
            for( let d of p.primary_site ){
                if( ! sites.includes(d) ){
                    sites.push(d);
                }
            }
        }
        
        this.meta.program = programs;
        this.meta.disease = diseases;
        this.meta.primary_site = sites;
    }

    async get_projects(){
        let r = await fetch('https://api.gdc.cancer.gov/projects?size=1000&sort=project_id:asc')
        let dat = await r.json();
        this.projects = dat.data.hits;
        this.pids = dat.data.hits.map( e => e.id );

        this.get_metadata_options();
    }
    
    async get_summary_by_project( proj_id ){
        let r = await fetch( `https://api.gdc.cancer.gov/projects/${proj_id}?expand=summary,summary.experimental_strategies,summary.data_categories` );
        let dat = await r.json();
        
        return dat.data;
    }
    
    async get_all_project_summary( proj_list ){
        let that = this;
        let dat = await Promise.all( proj_list.map( proj_id => that.get_summary_by_project( proj_id ) ) );
        let counts = {};
        let info = {};
        for (let el of dat ){
            info[ el.project_id ] = el.summary;
            counts[ el.project_id ] = { "cases": el.summary.case_count, "files": el.summary.file_count };
        }
        this.projects_summary = info;
        this.processed_counts = counts;
    }

    async _get_cases_demography(p){
        // cases.demographic.ethnicity, cases.demographic.gender, cases.demographic.race
        let filters = { "op":"and",
                "content":[
                    {
                        "op":"in",
                        "content":{
                            "field":"project.project_id",
                                "value": [ p ]
                        }
                    }
                ]
            }
        let body = JSON.stringify( { "filters": filters, "facets": "demographic.gender,demographic.ethnicity,demographic.race", "size": 0 } );
        let r = await fetch("https://api.gdc.cancer.gov/v0/cases", { "method": "POST", "headers": { "Content-Type": "application/json" }, "body": body });
        let dat = await r.json();

        return dat.data.aggregations;
    }

    async _get_cases_count_by_data_category(p){
        // cases.demographic.ethnicity, cases.demographic.gender, cases.demographic.race
        let filters = { "op":"and",
                "content":[
                    {
                        "op":"in",
                        "content":{
                            "field":"project.project_id",
                                "value": [ p ]
                        }
                    }
                ]
            }
        let body = JSON.stringify( { "filters": filters, "facets": "files.data_category", "size": 0 } );
        let r = await fetch("https://api.gdc.cancer.gov/v0/cases", { "method": "POST", "headers": { "Content-Type": "application/json" }, "body": body });
        let dat = await r.json();

        return dat.data.aggregations;
    }

    async _get_cases_count_by_exp_strategy(p){
        // cases.demographic.ethnicity, cases.demographic.gender, cases.demographic.race
        let filters = { "op":"and",
                "content":[
                    {
                        "op":"in",
                        "content":{
                            "field":"project.project_id",
                                "value": [ p ]
                        }
                    }
                ]
            }
        let body = JSON.stringify( { "filters": filters, "facets": "files.experimental_strategy", "size": 0 } );
        let r = await fetch("https://api.gdc.cancer.gov/v0/cases", { "method": "POST", "headers": { "Content-Type": "application/json" }, "body": body });
        let dat = await r.json();

        return dat.data.aggregations;
    }

    async prepare_coverage_data_dataCategory_expStrategy(p, datcat){
        let that = this;
        let promises = [ this._get_cases_count_by_data_category(p), this._get_cases_count_by_exp_strategy(p) ];
        let dat = await Promise.all( promises );

        let total_cases = this.processed_counts[p].cases;
        let coverage = { 
            "data_category": { 
                "x": dat[0]["files.data_category"].buckets.map( e => e.key ),
                "absolute": dat[0]["files.data_category"].buckets.map( e => e.doc_count ), 
                "relative": dat[0]["files.data_category"].buckets.map( e => ( e.doc_count / total_cases )*100 )
            }, 
            "experimental_strategy": { 
                "x": dat[1]["files.experimental_strategy"].buckets.map( e => e.key ),
                "absolute": dat[1]["files.experimental_strategy"].buckets.map( e => e.doc_count ), 
                "relative": dat[1]["files.experimental_strategy"].buckets.map( e => ( e.doc_count / total_cases )*100 )
            }
        };

        return coverage;
    }
    
    async get_case_files_by_data_category(p, datcat){
        let format = this.formats_datCategory[datcat];
        
        let filters = { "op":"and",
                "content":[
                    {
                        "op":"in",
                        "content":{
                            "field":"cases.project.project_id",
                                "value": [ p ]
                        }
                    },
                    {
                        "op":"=",
                        "content":{
                            "field":"data_category",
                                "value": datcat
                        }
                    },
                    {
                        "op":"=",
                        "content":{
                            "field":"access",
                                "value": "open"
                        }
                    },
                    {
                        "op":"in",
                        "content":{
                            "field":"data_format",
                                "value": format
                        }
                    }
                ]
            }
        filters = encodeURI( JSON.stringify( filters ) );
        
        let url = `https://api.gdc.cancer.gov/files?filters=${filters}&fields=file_id,file_size,platform,cases.project.project_id,cases.submitter_id,cases.case_id,cases.samples.tumor_descriptor,cases.samples.tissue_type,cases.demographic.ethnicity,cases.demographic.gender,cases.demographic.race,cases.demographic.year_of_birth,cases.diagnoses.vital_status,cases.diagnoses.days_to_last_follow_up,cases.diagnoses.age_at_diagnosis,cases.diagnoses.classification_of_tumor,cases.diagnoses.days_to_recurrence,cases.diagnoses.tumor_stage&size=1000`
        let r = await fetch( url );
        let dat = await r.json();

        let before = dat.data.hits.length;
        let that = this;
        dat = dat.data.hits.filter( e => e.file_size <= that.file_size_limit );
        let after = dat.length;
        let delta = after - before;
        if( delta > 0 ){
        console.log( `INFORMATION: ${delta} files were filtered out because their sizes are up to to the file size limit configured (${this.file_size_limit})` )
        }

        // stages = dat.data.hits.map( e => e.cases[0].samples[0]['tumor_descriptor'] )

        return dat;
    }
    
    async _get_file_by_uuid(uuid){
        let url = `https://api.gdc.cancer.gov/data/${uuid}`;
        let r = await fetch( url );
        let dat = await r.text();
        let lines = dat.split('\n').slice(0, -1);

        return lines;
        // values = dat.split('\n').map( e => { try{ return parseFloat( e.split('\t')[1] ) } catch { return 0 } } )
    }

    async get_clinical_stratification_survival(project){
        let url = `${location.href}/data_processed/${project}/clinical/data_cases.json`;
        let r = await fetch( url );
        let dat_cases = await r.json();

        url = `${location.href}/data_processed/${project}/clinical/survival_probs.json`;
        r = await fetch( url );
        let dat_survival = await r.json();

        let result = { "cases": dat_cases, "survival": dat_survival };

        return result;
    }
    
    async retrieve_process_methylation_files( uuids ){
        let that = this;
        this.metexp_features = {};
        let promises = uuids.map( id => that._process_methylation_file(id) );
        let dat = await Promise.all( promises );

        return dat;
    }

    async _process_methylation_file(uuid){
        let that = this;
        let lines = await this._get_file_by_uuid(uuid);
        let annotation = {};
        let cgids = Object.keys(this.methylation_cgids_ann);
        lines = lines.map( e => e.split('\t') ).filter( e => cgids.includes(e[0]) );
        dt = this.metexp_features;
        lines.forEach( (l) => {
            let keys = Object.keys(dt);
            if( ! keys.include(l[0]) ){
                dt[ l[0] ] = {};
            }
            dt[ l[0] ][uuid] = parseFloat(l[1]);
        } );
        this.metexp_features = dt;


        /*
        let dt = {};
        lines.forEach( (l) => {
            let ann = that.methylation_cgids_ann[l[0]]
            ann["value"] = parseFloat(l[1]);
            dt[l[0]] = ann;
        } );
        
        return dt;
        */
        return 1
    }

    async get_methylation_mapping(){
        let url = `${location.href}/data_processed/all_mapp.json`;
        let r = await fetch( url );
        let dat = await r.json();

        this.methylation_cgids_ann = dat;

        return dat;
    }
    
    async get_case_data_survival(p){
        
        let filters = { "op":"and",
                "content":[
                    {
                        "op":"in",
                        "content":{
                            "field":"project.project_id",
                                "value": [ p ]
                        }
                    }
                ]
            }
        filters = encodeURI( JSON.stringify( filters ) );
        
        let url = `https://api.gdc.cancer.gov/v0/cases?filters=${filters}&fields=case_id,diagnoses.age_at_diagnosis,diagnoses.classification_of_tumor,diagnoses.days_to_death,diagnoses.days_to_last_follow_up,diagnoses.classification_of_tumor,files.cases.diagnoses.vital_status,demographic.gender,file.cases.diagnoses.tumor_stage&size=1000`
        let r = await fetch( url );
        let dat = await r.json();

        // stages = dat.data.hits.map( e => e.cases[0].samples[0]['tumor_descriptor'] )

        return dat.data;
    }

    async map_fileid_to_caseid(project, datcat){
        datcat = datcat.replaceAll(" ", "-");
        let url = `${location.href}/data_processed/${project}/${datcat}/files_metadata.tsv`;
        let r = await fetch( url );
        let dat = await r.text();
        let txt = dat.split("\n").slice(1, -1).map( l => l.split('\t') );
        let mapp = {};
        txt.forEach( l => { mapp[ l.slice(-2)[0] ] = l[0]; } );

        return mapp;
    }

    async map_caseid_to_demovars(project){
        let datcat = "cases";
        datcat = datcat.replaceAll(" ", "-");
        let url = `${location.href}/data_processed/${project}/${datcat}/cases_metadata.json`;
        let r = await fetch( url );
        let mapp = await r.json();

        return mapp;
    }

    async get_mutation_info_metrics_stratified(p){
        let datcat = "simple nucleotide variation";
        datcat = datcat.replaceAll(" ", "-");
        let cols_stratification = ['race','gender', 'ethnicity'];
        let cols_info = ["impact", "consequence", "clin_sig", "variant_class", "hugo_symbol"];
        let mp_name_info = { "impact": "Impact", "consequence": "Consequence", "clin_sig": "Pathogenicity", "hugo_symbol": "Gene" };
        // map nan to 'Not reported'
        
        let mapf = await this.map_fileid_to_caseid(p, datcat);
        let mapc = await this.map_caseid_to_demovars(p);
        let url = `${location.href}/data_processed/${p}/${datcat}/general_snprs_info.json`;
        let r = await fetch( url );
        let dat = await r.json();

        let rdat_mut = {}; // subgroup -> mutation -> count
        // section 1 plot in fact will be rdat_mut['by_gender'] = { 'male': { 'mapk1_v200f': 3 } } | by_gender or by_colStratification (second select in page) -> in_common or exclusive (first section) -> subgroup , value of col_stratification (male or female) -> mutation -> count
        // Provide a full table ?

        let rdat_mdetails = {}; // by_gender or by_colStratification (second select in page) -> col_info (impact, consequence, etc) (different plot area) -> subgroup -> values col_info (impacto LOW, HIGH, etc) -> count
        // This one can be grouped, but the col_info values that are missing must be filled with zero

        // subgroup -> mutation
        let all_muts = {};
        let all_values_mut_subgroup = {}; // to calculate later in_common and exclusive
        let all_values_details = {};
        for(let c of cols_stratification){
            all_values_mut_subgroup[`by_${c}`] = {};
            rdat_mdetails[`by_${c}`] = {};
        }

        // Getting all possible values first
        let uuids = Object.keys(dat);
        for( let f of uuids){
            let case_id = mapf[f];
            let infosub = mapc[case_id];

            let mutations = Object.keys(dat[f]);
            for(let m of mutations){
                let infom = dat[f][m];
                if( ! Object.keys(all_muts).includes(m) ){
                    all_muts[m] = 0;
                }
                all_muts[m] += 1;

                for ( let ci of cols_info){
                    if( ! Object.keys(all_values_details).includes(ci) ){
                        all_values_details[ci] = new Set();
                    }
                    all_values_details[ci].add( infom[ci] );
                }

                
                for(let c of cols_stratification){
                    let sub = infosub[c];
                    if( ! Object.keys(all_values_mut_subgroup[`by_${c}`]).includes(sub) ){
                        all_values_mut_subgroup[`by_${c}`][sub] = {};
                        rdat_mdetails[`by_${c}`][sub] = {}
                    }

                    if( ! Object.keys(all_values_mut_subgroup[`by_${c}`][sub]).includes(m) ){
                        all_values_mut_subgroup[`by_${c}`][sub][m] = 0;
                    }

                    all_values_mut_subgroup[`by_${c}`][sub][m] += 1;

                    for ( let ci of cols_info){
                        if( ! Object.keys( rdat_mdetails[`by_${c}`][sub] ).includes(ci) ){
                            rdat_mdetails[`by_${c}`][sub][ci] = {};
                        }
                        if( ! Object.keys( rdat_mdetails[`by_${c}`][sub][ci] ).includes( infom[ci] ) ){
                            rdat_mdetails[`by_${c}`][sub][ci][ infom[ci] ] = 0;
                        }
                        rdat_mdetails[`by_${c}`][sub][ci][ infom[ci] ] += 1;
                    }
                }
            }
        }

        // Composing the final data
        let in_common = {};
        let exclusive = {};
        for(let c of cols_stratification){
            in_common[`by_${c}`] = _get_intersection_values( all_values_mut_subgroup[`by_${c}`], "dict");

            exclusive[`by_${c}`] = {};
            let groups = Object.keys( all_values_mut_subgroup[`by_${c}`] );
            for( let sub of groups ){
                exclusive[`by_${c}`][sub] = _get_exclusive_values(sub, all_values_mut_subgroup[`by_${c}`], "dict");

                for( let ci of cols_info ){
                    for( let v of all_values_details[ci] ){
                        if( ! Object.keys( rdat_mdetails[`by_${c}`][sub][ci] ).includes(v) ){
                            rdat_mdetails[`by_${c}`][sub][ci][ v ] = 0;
                        }
                    }
                }
            }
        }
        console.log("all_muts", all_muts);
        console.log("in_common", in_common);
        console.log("exclusive", exclusive);
        console.log("all_values_mut_subgroup", all_values_mut_subgroup);
        console.log("rdat_mdetails", rdat_mdetails);
        console.log("all_values_details", all_values_details);

        return [all_muts, in_common, exclusive, all_values_mut_subgroup, rdat_mdetails, all_values_details];
    }

    async load_cases_mutation_annotations(project){
        // Cases
        let url = `${location.href}/data_processed/${project}/clinical/data_cases.json`;
        let r = await fetch( url );
        let dat_cases = await r.json();

        // Mutation information
        let datcat = "simple nucleotide variation";
        datcat = datcat.replaceAll(" ", "-");
        let cols_stratification = ['race','gender', 'ethnicity'];

        let dm = {};
        for(let c of cols_stratification){
            dm[`by_${c}`] = {};
        
            let url = `${location.href}/data_processed/${project}/${datcat}/by_${c}_table_summary.tsv`;
            let r = await fetch( url );
            let txt = await r.text();
            let dat = _format_table(txt, ["count"]);
            for(let el of dat){
                let sub = el["subgroup"];
                let ci = el["feature"];
                let v = el["feature_value"];
                let cnt = el["count"];
                
                if( ! Object.keys( dm[`by_${c}`] ).includes(sub) ){
                    dm[`by_${c}`][sub] = {}
                }
                if( ! Object.keys( dm[`by_${c}`][sub] ).includes(ci) ){
                    dm[`by_${c}`][sub][ci] = {}
                }
                
                dm[`by_${c}`][sub][ci][v] = cnt;
            }
        }

        let result = { "cases": dat_cases, "ann_mutations": dm }
        return result;
    }

    
}
