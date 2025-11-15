class CoverageSummary extends HTMLElement {
  constructor() {
    super();
  }

  connectedCallback() {
    this.innerHTML = `
      <p id = "notice" > Loading... </p>

      <section id="coverage_analysis" class="mt-3">
            <p>
                Overview about the GDC TCGA projects available and the cases coverage by omic data.
            </p>
            
            <div id="project_filter" style="display: none;" > 
                <h4> Filter the projects you want to explore (max. 40, out of <span id="total_projects">  </span> ) </h4>

                <div class="row g-2 align-items-center" id="filters_area" >


                </div>

                <div class="mt-3 col-12" >
                    <button type="button" class="ml-3 btn btn-primary" id="go_filter" onclick="get_filtered_projects()" > Filter projects </button>

                    <button type="button" class="ml-3 btn btn-primary" id="go_filter" onclick="update_main_plots()" > Reload plots </button>
                </div>

                <div class="mt-3 col-12" >
                    <h5> Selected projects: </h5>
                    <div id="projects_set" > </div>
                </div>
            </div>

            <hr />

            <div class="row justify-content-center mt-3"  >
                <div class="col-12" id="instructions" style="display: none;" >
                    <p>
                        Choose a project (clicking on the blue bar) to load details about the information available for the cases.
                    </p>

                    <div  class="mt-3 row justify-content-center"  >
                        <div id="plot_cases_files_count" class="col-12 " style="display: none;" >  </div>
                    </div>
                </div>

                <div class="col-12" id="analysis_current" style="display: none;" >
                    <h4 class="mt-3" id="selected_proj" >  </h4>
                    <p>
                        <span style = "font-weight: bold;" > Name: </span> 
                        <span id = "proj_name" >  </span> 
                    </p>

                    <!-- plot demographic information of cases -->
                    <div class="col-12" id="cases_demography" style="display: none;" >
                        <p> The plots below show the distribution of cases across demographic information </p>
                        <div id="plot_cases_demography" class="mt-3 row justify-content-start"  >  </div>
                    </div>

                    <!-- plot coverage -->
                    <div class="col-12" id="cases_coverage" style="display: none;" >
                        <p> The plots below show the cases coverage for the technical information regarding the project files repository </p>
                        
                        <div  class="mt-3 row justify-content-start"  >
                            <div id="plot_cases_coverage_datCategory" class="mt-3 col-6" >  </div>

                            <div id="plot_cases_coverage_expStrategy" class="mt-3 col-6" >  </div>
                        </div>
                    </div>
                </div>

            </div>
        </section>
        
    `;
  }
}
customElements.define('coverage-component', CoverageSummary);

/*
Page organization:
sec 1 - grouped bar plot, x axis are the projects and the bar colors are the cases and files count
    when click in file count bar of a project show pie plot with those that are public or restrict.
sec 2 - a text field for the person filter the projects based on the disease or primary site key word
    sec 2 - there is a select to filter the project.
    plot showing the coverage in % of data categories cases over total cases of proj
    upon click on data cat. bar it shows in a pie plot the coverage in exp strategy available
        document.getElementById("myDiv").on('plotly_click', function(data){
            console.log(data.points[0].x);
        });
sec 3 - demographic arrangement of cases. given a project (select) there is also a select to filter the view by race, gener, ethnicity, age range (baby, child, teen, adult, elderly)
*/

let obj_cov = {};

let config = { "responsive": true };
let pie_layout = { "width": 500, "height": 500 };

// Section 0 - projects filtering
// obj_cov.projects.filter( e => e.disease_type.map( e => e.toLowerCase() ).filter( it => it.indexOf('epith')!=-1 ) )

function render_filter_projs_area(){
    project_filter.style.display = 'none';

    let domid_container = "filters_area";

    let selected = "TCGA";
    let domid_target = "program";
    let label = _capitalize(domid_target);
    let options = obj_cov.meta.program;
    fill_select( label, options, domid_target, domid_container, selected, 'Any');

    selected = "";
    domid_target = "disease";
    label = _capitalize(domid_target);
    options = obj_cov.meta.disease;
    options.unshift('Any');
    fill_select( label, options, domid_target, domid_container, selected, 'Any');

    selected = "";
    domid_target = "primary_site";
    label = _capitalize(domid_target);
    options = obj_cov.meta.primary_site;
    options.unshift('Any');
    fill_select( label, options, domid_target, domid_container, selected, 'Any');

    total_projects.innerHTML = obj_cov.projects.length;

    project_filter.style.display = '';
}

function update_main_plots(){
    analysis_current.style.display = 'none';
    instructions.style.display = 'none';

    let p = select_program.value;
    let d = select_disease.value;
    let s = select_primary_site.value;

    let fp = new Set( obj_cov.projects.filter( x => (x.project_id.split('-')[0] == p) || p == 'Any' ) );
    let fd = new Set( obj_cov.projects.filter( x => (x.disease_type.includes(d) ) || d == 'Any' ) );
    let fs = new Set( obj_cov.projects.filter( x => (x.primary_site.includes(s) ) || s == 'Any' ) );
    
    let inter = intersect(fp, fd, fs);
    obj_cov.filtered_projects = Array.from(inter).map( x => x.project_id );

    generate_tags(obj_cov.filtered_projects, "projects_set");

    prepare_plot_overview_projects();
    
    instructions.style.display = '';

}

function _render_pieplot_demography(){
    let p = obj_cov.current_project.project_id;

    obj_cov._get_cases_demography(p).then( (dat) => {
        let keys = Object.keys(dat);
        let div_ids = [];
        let htmls = "";
        keys.forEach( (it) => {
            let _id = it.split('.')[1];

            div_ids.push(`pie_${_id}`);
            htmls += `
                <div  id = "pie_${_id}" class="col-4" >

                </div>
            `;
        });
        document.getElementById("plot_cases_demography").innerHTML = htmls;

        keys.forEach( (it) => {
            let _id = it.split('.')[1];

            let itlay = pie_layout;
            itlay["title"] = { "text": `By ${ _id }` };
            let labels = dat[it].buckets.map( e => e.key );
            let values = dat[it].buckets.map( e => e.doc_count );
            let pldata = [{ "values": values, "labels": labels, "type": "pie" }];
            Plotly.newPlot( `pie_${_id}`, pldata, itlay, config);
        });
    
        if( keys.length > 0 ){
            document.getElementById("cases_demography").style.display = '';
        }
    } );
}

function setup_current_project_analysis(p){
    analysis_current.style.display = 'none';

    obj_cov.current_project = obj_cov.projects.filter( x => x.project_id == p )[0];
    document.getElementById('proj_name').innerHTML = obj_cov.current_project.name;

    selected_proj.innerHTML = `Selected project: ${p}`;

    _render_pieplot_demography();
    prepare_plots_coverage();

    analysis_current.style.display = '';  
}

// Section 1
function prepare_plot_overview_projects(){

    let layout = { barmode: 'group', title: { text: "Cases and Files count per project" }, xaxis: { title: { text: 'Projects' } }, yaxis: { title: { text: 'Cases Count' } }, yaxis2: { title: 'Files Count', overlaying: 'y', side: 'right' } };
    
    let filtered_pids = obj_cov.filtered_projects; // obj_cov.pids.filter( e => e.toLowerCase().indexOf('tcga')!=-1 );
    
    let pl = [];
    let y = filtered_pids.map( e => obj_cov.processed_counts[e].cases );
    pl.push( { "name": "Cases", "x": filtered_pids, "y": y, "type": "bar" } );

    y = filtered_pids.map( e => obj_cov.processed_counts[e].files );
    pl.push( { "name": "Files", "x": filtered_pids, "y": y, "type": "scatter", "mode": 'lines+markers', "yaxis": 'y2' } );
    
    let sec1plot = document.getElementById('plot_cases_files_count');
    Plotly.newPlot( "plot_cases_files_count", pl, layout, config );

    sec1plot.on('plotly_click', function(data){
        let p = data.points[0].x;
        console.log(p);
        setup_current_project_analysis(p);
        
    });
    sec1plot.style.display = '';
}
    
// Section 2
function _render_plot_cases_coverage_datCategory(info){
    let layout = { 
        barmode: 'group', 
        title: { text: "Cases coverage per Data Category" }, 
        xaxis: { title: { text: 'Data Category' } }, 
        yaxis: { title: { text: 'Abs. Cases Count' } }, 
        yaxis2: { title: 'Cases Coverage (%)', overlaying: 'y', side: 'right', range: [0, 100] }
    };
    
    let pl = [];
    let y = info.absolute;
    pl.push( { "name": "Absolute Cases count", "x": info.x, "y": y, "type": "bar" } );

    y = info.relative;
    pl.push( { "name": "Cases Coverage (%)", "x": info.x, "y": y, "type": "scatter", "mode": 'lines+markers', "yaxis": 'y2' } );
    
    let sec2plot1 = document.getElementById('plot_cases_coverage_datCategory');
    Plotly.newPlot( "plot_cases_coverage_datCategory", pl, layout, config );

    sec2plot1.style.display = '';
}
function _render_plot_cases_coverage_expStrategy(info){
    let layout = { 
        barmode: 'group', 
        title: { text: "Cases coverage per Experimental Strategy" }, 
        xaxis: { title: { text: 'Data Category' } }, 
        yaxis: { title: { text: 'Abs. Cases Count' } }, 
        yaxis2: { title: 'Cases Coverage (%)', overlaying: 'y', side: 'right', range: [0, 100] } 
    };
    
    let pl = [];
    let y = info.absolute;
    pl.push( { "name": "Absolute Cases count", "x": info.x, "y": y, "type": "bar" } );

    y = info.relative;
    pl.push( { "name": "Cases Coverage (%)", "x": info.x, "y": y, "type": "scatter", "mode": 'lines+markers', "yaxis": 'y2' } );
    
    let sec2plot2 = document.getElementById('plot_cases_coverage_expStrategy');
    Plotly.newPlot( "plot_cases_coverage_expStrategy", pl, layout, config );

    sec2plot2.style.display = '';
}

function prepare_plots_coverage(){
    cases_coverage.style.display = 'none';

    let p = obj_cov.current_project.project_id;

    obj_cov.prepare_coverage_data_dataCategory_expStrategy(p).then( (dat) => {
        _render_plot_cases_coverage_datCategory( dat["data_category"] );
        _render_plot_cases_coverage_expStrategy( dat["experimental_strategy"] );

        cases_coverage.style.display = '';
    } );
}

function render_coverage_by_project(){
  document.getElementById('project_select').innerHTML = 'Loading ...';
  
  let options = "";
  obj_cov.pids.forEach( (v) => { options += `<option value="${v}" > ${v} </option>`; } );
  
  let html = `
  <label class="form-label" > Project:</label>
  <select class="form-control" id="sel_project">
      ${options}
  </select>
  `;
  document.getElementById('project_select').innerHTML = html;
  
  render_plots_cpu_alloc();
}

let init_case_coverage = async () => {
    document.getElementById('notice').innerHTML = 'Loading...';
      
    obj_cov = new GdcExplorerLib( );
    await obj_cov.get_projects();
    await obj_cov.get_all_project_summary( obj_cov.pids );
    
    render_filter_projs_area();
    update_main_plots();

    document.getElementById('notice').innerHTML = '';
}

init_case_coverage().then( async v => { console.log("Initialized!"); } );


