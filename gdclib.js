class GdcExplorerLib {

    constructor(){
        this.pids = [];
        this.projects = [];
    	this.projects_summary = {};
    	this.processed_counts = {};
    	
    }
    // Coverage is given by the case_couunt in each entry of exp strategy or data category divided by the general case_count of the project
    
    async get_projects(){
        let r = await fetch('https://api.gdc.cancer.gov/projects?size=1000&sort=project_id:asc')
        let dat = await r.json();
        this.projects = dat.data.hits;
        this.pids = dat.data.hits.map( e => e.id );
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
    
}
