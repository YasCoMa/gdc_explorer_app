class ExpAiDisparitySummary extends HTMLElement {
  constructor() {
    super();
  }

  connectedCallback() {
    this.innerHTML = `
      <p id = "notice_ai1" > Loading... </p>

      <section id="coverage_analysis" class="mt-3">
            <p>
                Exploratory data analysis and feature extraction with group stratification, acccording to a project and a target data category
            </p>
            
            <div id="project_filter_ai1" style="display: none;" > 
                <h4> Filter the project you want to explore </h4>

                <div class="row g-2 align-items-center" id="filters_area_ai1" >


                </div>

                <div class="mt-3 col-12" >
                    <button type="button" class="ml-3 btn btn-primary" id="go_filter" onclick="analyze_eda()" > Analyze </button>
                </div>
            </div>

            <hr />

            <div class="row justify-content-center mt-3"  >
                <div class="col-12" id="instructions_ai1" style="display: none;" >
                    <p>
                        Choose a project (clicking on the blue bar) to load details about the information available for the cases.
                    </p>
                </div>

                <div class="col-12" id="analysis_current_ai1" style="display: none;" >
                    
                </div>

            </div>
        </section>
        
    `;
  }
}
customElements.define('expai1-component', ExpAiDisparitySummary);

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

let obj_ai1 = {};

let config_ai1 = { "responsive": true };
let pie_layout_ai1 = { "width": 500, "height": 500 };

// Section 0 - projects filtering
// obj_ai1.projects.filter( e => e.disease_type.map( e => e.toLowerCase() ).filter( it => it.indexOf('epith')!=-1 ) )

function _fill_data_category(){
    let p = select_project_ai1.value;
    let allowed = new Set( Object.keys(obj_ai1.formats_datCategory) )

    let domid_container = "filters_area_ai1";

    let domid_target = "data_category_ai1";
    let label = "Data Category";
    let options = new Set( obj_ai1.projects_summary[p].data_categories.map( e => e.data_category.toLowerCase() ) );
    options = options.intersection(allowed); // Filter data categories that have open parsable (not platform-specific raw output as idat) data
    options = Array.from( options );
    let selected = options[0];
    fill_select( label, options, domid_target, domid_container, selected, null);
}

function render_filter_projs_area_ai1(){
    project_filter_ai1.style.display = 'none';

    let domid_container = "filters_area_ai1";

    let domid_target = "project_ai1";
    let label = "Project";
    let options = Object.keys(obj_ai1.projects_summary);
    
    let selected = options[0];
    let onchange = "_fill_data_category()"
    fill_select( label, options, domid_target, domid_container, selected, onchange);

    _fill_data_category();

    project_filter_ai1.style.display = '';
}

let init_case_expai1 = async () => {
    document.getElementById('notice_ai1').innerHTML = 'Loading...';
      
    obj_ai1 = new GdcExplorerLib( );
    await obj_ai1.get_projects();
    await obj_ai1.get_all_project_summary( obj_ai1.pids );
    
    render_filter_projs_area_ai1();

    document.getElementById('notice_ai1').innerHTML = '';
}

init_case_expai1().then( async v => { console.log("Initialized!"); } );


