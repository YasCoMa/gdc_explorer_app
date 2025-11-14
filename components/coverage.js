class CoverageSummary extends HTMLElement {
  constructor() {
    super();
  }

  connectedCallback() {
    this.innerHTML = `
      <section id="total_cpus_analysis" class="mt-3">
          <p>
             Overview about the GDC projects available and the cases coverage by omic data.
          </p>
          
          <p>
             Choose below if you want to see the users usage in each partition or in each node. You can also specify a specific one, if you want.
          </p>
          
          <div class="row g-2 align-items-center">
            <div class="col-auto">
                <label class="form-label" > View options:</label>
                <select class="form-control" id="view_cpu_alloc" onchange="update_context_cpu_alloc()" >
                  <option value="partition" selected > Per Partition </option>
                  <option value="node" > Per Node </option>
                </select>
            </div>
            
            <div class="col-auto" id="specific_filter_cpu_alloc" onchange="render_plots_cpu_alloc()" >
            </div>
            
            <div class="col-auto" >
                <button type="button" class="ml-3 btn btn-primary" id="go_cpu_alloc" onclick="render_plots_cpu_alloc()" >Reload plots</button>
            </div>
          </div>
          
          <div id="analysis_part_cpu_alloc" class="mt-3" style="display: none;" >  </div>
          
        </section>
        
        <hr />
        
    `;
  }
}
customElements.define('coverage-component', CoverageSummary);

/*
Page organization:
sec 1 - grouped bar plot, x axis are the projects and the bar colors are the cases and files count
sec 2 - a text field for the person filter the projects based on the disease or primary site key word
    sec 2 - there is a select to filter the project.
    plot showing the coverage in % of data categories cases over total cases of proj
    upon click on data cat. bar it shows in a pie plot the coverage in exp strategy available
        document.getElementById("myDiv").on('plotly_click', function(data){
            console.log(data.points[0].x);
        });
sec 3 - demographic arrangement of cases. given a project (select) there is also a select to filter the view by race, gener, ethnicity, age range (baby, child, teen, adult, elderly)
*/

let obj_cpu_alloc = {};
                  
function update_context_cpu_alloc(){
  let view = document.getElementById('view_cpu_alloc').value;
  
  let name = _capitalize(view);
  let values = obj_cpu_alloc[view+'s'];
  let options = "";
  values.forEach( (v) => { options += `<option value="${v}" > ${v} </option>`; } );
  
  let html = `
  <label class="form-label" > ${name}:</label>
  <select class="form-control" id="view_filter_cpu_alloc">
      <option value="all" > All </option>
      ${options}
  </select>
  `;
  document.getElementById('specific_filter_cpu_alloc').innerHTML = html;
  
  render_plots_cpu_alloc();
}

function render_plots_cpu_alloc(){
  document.getElementById('analysis_part_cpu_alloc').innerHTML = "";
  go_cpu_alloc.disabled = true;
  go_cpu_alloc.innerHTML = "Loading...";
      
  let new_div_els = "";
  let new_div_ids = [];
  let layout = { barmode: 'group', xaxis: { title: { text: 'Weeks of last 2 months' } }, yaxis: { title: { text: 'Average CPUs allocated for jobs' } } };
  
  let view = document.getElementById('view_cpu_alloc').value;
  let filter = document.getElementById('view_filter_cpu_alloc').value;
  let name = _capitalize(view);
  
  let dat = obj_cpu_alloc.prepare_plot_data_alloc_cpu (view, filter);
  
  if( Object.keys( dat ).length > 0 ){
    let j = 1;
    for( let v of Object.keys( dat ) ){
      let _id = `cpualloc_${view}_plot_${j}`;
      new_div_ids.push(_id);
      new_div_els += `<div id="${_id}" > </div>`;
      j +=1 ;
    }
    document.getElementById('analysis_part_cpu_alloc').innerHTML = new_div_els;
    
    j = 0;
    for( let v of Object.keys(dat) ){
      let pldat = dat[v];
      let itlay = layout;
      itlay["title"] = { text: `${name} - ${v}` };
      Plotly.newPlot( new_div_ids[j], pldat, itlay );
      j += 1;
    }
  }
  else{
    document.getElementById('analysis_part_cpu_alloc').innerHTML = "<p> No data available at this moment! </p>";
  }
  
  document.getElementById('analysis_part_cpu_alloc').style.display = '';
  go_cpu_alloc.disabled = false;
  go_cpu_alloc.innerHTML = "Reload plots";
  
}

let init_case_cpu_alloc = async () => {
      go_cpu_alloc.disabled = true;
      go_cpu_alloc.innerHTML = "Loading...";
      
      obj_cpu_alloc = new IrbclusterJobMonitor( );
      await obj_cpu_alloc.load_data_alloc('cpu');
      update_context_cpu_alloc();
      render_plots_cpu_alloc();
      
      go_cpu_alloc.disabled = false;
      go_cpu_alloc.innerHTML = "Reload plots";
}

init_case_cpu_alloc().then( async v => { console.log("Initialized!"); } );


