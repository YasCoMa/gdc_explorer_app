class CoverageSummary extends HTMLElement {
  constructor() {
    super();
  }

  connectedCallback() {
    this.innerHTML = `
      <section id="total_cpus_analysis" class="mt-3">
          <p>
             Overview about the GDC TCGA projects available and the cases coverage by omic data.
          </p>
          
          <!--
          <p>
             Choose a project to load details about the information available for the cases.
          </p>
            
            <div class="col-auto" id="project_select" onchange="render_plots_cpu_alloc()" >
            </div>
            
            <div class="col-auto" >
                <button type="button" class="ml-3 btn btn-primary" id="go_cpu_alloc" onclick="render_plots_cpu_alloc()" >Reload plots</button>
            </div>
          </div>
          -->
          
          <div id="plot_cases_files_count" class="mt-3" style="display: none;" >  </div>
          
        </section>
        
        <hr />
        
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
         
// Section 1
function prepare_plot_overview_projects(){
    let layout = { barmode: 'group', title: { text: "Cases and Files count per project" }, xaxis: { title: { text: 'Projects' } }, yaxis: { title: { text: 'Cases Count' } }, yaxis2: { title: 'Files Count', overlaying: 'y', side: 'right' } };
    
    let tcga_pids = obj_cov.pids.filter( e => e.toLowerCase().indexOf('tcga')!=-1 );
    
    let pl = [];
    let y = tcga_pids.map( e => obj_cov.processed_counts[e].cases );
    pl.push( { "name": "Cases", "x": tcga_pids, "y": y, "type": "bar" } );
    y = tcga_pids.map( e => obj_cov.processed_counts[e].files );
    pl.push( { "name": "Files", "x": tcga_pids, "y": y, "type": "scatter", yaxis: 'y2' } );
    
    Plotly.newPlot( "plot_cases_files_count", pl, layout );
    document.getElementById('plot_cases_files_count').style.display = '';
}
    
// Section 2              
function render_project_selection(){
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
      
      obj_cov = new GdcExplorerLib( );
      await obj_cov.get_projects();
      await obj_cov.get_all_project_summary( obj_cov.pids );
      
      prepare_plot_overview_projects();
}

init_case_coverage().then( async v => { console.log("Initialized!"); } );


