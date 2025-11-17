class GdcExplorerLib {

    constructor(){
        this.meta = { 'program': [], 'disease': [], 'primary_site': [] };
        this.pids = [];
        this.projects = [];
    	this.projects_summary = {};
    	this.processed_counts = {};
    	
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

    async prepare_coverage_data_dataCategory_expStrategy(p){
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
        }

        return coverage;
    }
    
}
