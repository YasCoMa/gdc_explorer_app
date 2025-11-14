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
          
            <div class="row justify-content-center">
                <div class="col-12">
                    <div id="plot_cases_files_count" class="mt-3 " style="display: none;" >  </div>
                </div>
                <div class="col-12">
                    <div id="plot_cases_demography" class="mt-3 " style="display: none;" >  </div>
                </div>
                <div class="col-6">
                    <div id="plot_files_access" class="mt-3 " style="display: none;" >  </div>
                </div>
            </div>
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

let config = { "responsive": true };
let pie_layout = { "width": 400, "height": 400 };


function render_pieplot_demography(p){
    obj_cov._get_cases_demography(p).then( (dat) => {
        let keys = Object.keys(dat);
        let div_ids = [];
        let htmls = "";
        keys.forEach( (it) => {
            let _id = it.split('.')[1];

            div_ids.push(`pie_${_id}`);
            htmls += `
                <div  id = "pie_${_id}" >

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
        document.getElementById("plot_cases_demography").style.display = '';
    } );
}
 
async function _get_files_access(p){
    let filters = { "op":"and",
            "content":[
                {
                    "op":"in",
                    "content":{
                        "field":"cases.project.project_id",
                            "value": [ p ]
                    }
                },
                /*,
                {
                    "op":"=",
                    "content":{
                        "field":"access",
                        "value":"open"
                    }
                }
                {
                    "op":"=",
                    "content":{
                        "field":"data_category",
                        "value":"Simple Nucleotide Variation"
                    }
                },
                {
                    "op":"=",
                    "content":{
                        "field":"data_type",
                        "value":"Masked Somatic Mutation"
                    }
                },
                {
                    "op":"=",
                    "content":{
                        "field":"data_format",
                        "value":"maf"
                    }
                },
                {
                    "op":"=",
                    "content":{
                        "field":"experimental_strategy",
                        "value":"WXS"
                    }
                }
                */
            ]
        }
    let body = JSON.stringify({ "filters": filters, "facets": "access", "size": 0 });
    let r = await fetch("https://api.gdc.cancer.gov/v0/files", { "method": "POST", "headers": { "Content-Type": "application/json" }, "body": body });
    let dat = await r.json();
    
    /*
    let query = encodeURI( JSON.stringify(filters) );
    //let url = "https://api.gdc.cancer.gov/files?filters="+query+"&format=tsv&fields=file_id,cases.project.project_id,cases.submitter_id,cases.case_id,cases.samples.tumor_descriptor,cases.samples.tissue_type,cases.demographic.ethnicity,cases.demographic.gender,cases.demographic.race,cases.demographic.year_of_birth,cases.diagnoses.vital_status,cases.diagnoses.days_to_last_follow_up,cases.diagnoses.age_at_diagnosis,cases.diagnoses.classification_of_tumor,cases.diagnoses.days_to_recurrence,cases.diagnoses.tumor_stage&size=1000";
    let url = "https://api.gdc.cancer.gov/files?filters="+query+"&format=tsv&fields=file_id&size=10000";
    let r = await fetch(url);
    let dat = await r.text();
    */
}

// Section 1
function prepare_plot_overview_projects(){

    let layout = { barmode: 'group', title: { text: "Cases and Files count per project" }, xaxis: { title: { text: 'Projects' } }, yaxis: { title: { text: 'Cases Count' } }, yaxis2: { title: 'Files Count', overlaying: 'y', side: 'right' } };
    
    let tcga_pids = obj_cov.pids.filter( e => e.toLowerCase().indexOf('tcga')!=-1 );
    
    let pl = [];
    let y = tcga_pids.map( e => obj_cov.processed_counts[e].cases );
    pl.push( { "name": "Cases", "x": tcga_pids, "y": y, "type": "bar" } );
    y = tcga_pids.map( e => obj_cov.processed_counts[e].files );
    pl.push( { "name": "Files", "x": tcga_pids, "y": y, "type": "scatter", "mode": 'lines+markers', "yaxis": 'y2' } );
    
    let sec1plot = document.getElementById('plot_cases_files_count');
    Plotly.newPlot( "plot_cases_files_count", pl, layout, config )

    sec1plot.on('plotly_click', function(data){
        let p = data.points[0].x;
        console.log(p);
        render_pieplot_demography(p);
    });
    sec1plot.style.display = '';
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
    document.getElementById('notice').innerHTML = 'Loading...';
      
    obj_cov = new GdcExplorerLib( );
    await obj_cov.get_projects();
    await obj_cov.get_all_project_summary( obj_cov.pids );
      
    prepare_plot_overview_projects();

    document.getElementById('notice').innerHTML = '';
}

init_case_coverage().then( async v => { console.log("Initialized!"); } );


