class GdcExplorerLib {

    constructor(){
        this.file_size_limit = 15000000;
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
        console.log(dat);

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
    
    async retrieve_process_methylation_files( uuids ){
        let that = this;
        let promises = uuids.map( id => that._process_methylation_file(id) );
        let dat = await Promise.all( promises );

        return dat;
    }

    async _process_methylation_file(uuid){
        let that = this;
        let lines = await this._get_file_by_uuid(uuid);
        let annotation = {};
        let cgids = Object.keys(this.methylation_cgids_ann);
        lines = lines.map( e => e.split('\t') ).filter( e => cgids.includes(e[0]) )
        let dt = {};
        lines.forEach( (l) => {
            let ann = that.methylation_cgids_ann[l[0]]
            ann["value"] = parseFloat(l[1]);
            dt[l[0]] = ann;
        } );
        return dt
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

    
}
