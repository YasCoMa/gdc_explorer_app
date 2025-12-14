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
                    <button type="button" class="ml-3 btn btn-primary" id="go_filter" onclick="setup_current_project_ai1()" > Analyze </button>
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
                    <h4 class="mt-3" id="selected_proj_ai1" >  </h4>
                    <p>
                        <span style = "font-weight: bold;" > Name: </span> 
                        <span id = "proj_name_ai1" >  </span> 
                    </p>

                    <p>
                        <span style = "font-weight: bold;" > Tumor tissue site: </span> 
                        <span id = "tissue_ai1" >  </span> 
                    </p>

                    <p>
                        <span style = "font-weight: bold;" > Histological type: </span> 
                        <span id = "hist_type_ai1" >  </span> 
                    </p>

                    <p>
                        <span style = "font-weight: bold;" > Subgroups count: </span> 
                        <span id = "samples_info_ai1" >  </span> 
                    </p>

                    <!-- plot stage per subgroup -->
                    <div class="col-12" id="cases_stage_ai1" style="display: none;" >
                        <p> The plots below show the distribution of pathological state in each subgroup </p>
                        <div id="plot_cases_hist_stage_ai1" class="mt-3 row justify-content-start"  >  </div>
                    </div>

                    <!-- plot drugs used per subgroup -->
                    <div class="col-12" id="cases_most_used_drugs_ai1" style="display: none;" >
                        <p> The plots below show the most frequently prescribed drugs in each subgroup </p>
                        <div id="plot_cases_prescribed_drugs_ai1" class="mt-3 row justify-content-start"  >  </div>
                    </div>

                    <!-- plot distribution continuous metric variables per subgroup (age at diagnosis, admin drugs and radiation periods, dosis_drugs, dosis_radiation, duration_treat_drug, duration_treat_radiation) / https://plotly.com/javascript/violin/ -->
                    <div class="col-12" id="cases_dist_ai1" style="display: none;" >
                        <p> The plots below show the distribution of the following continuous variables in each subgroup: 
                            <ul>
                                <li>Age at diagnosis</li>
                                <li>Qty. followup events</li>
                                <li>Qty. drugs administration</li>
                                <li>Qty. radiation prescriptions</li>
                                <li>Duration of treatment by drugs</li>
                                <li>Duration of treatment by radiation</li>
                                <li>Dosis used in treatment by drugs</li>
                                <li>Dosis used in treatment by radiation</li>
                            </ul>
                        </p>
                        <div id="plot_cases_dist_ai1" class="mt-3 row justify-content-start"  >  </div>
                    </div>

                    <!-- plot survival KM per subgroup / https://plotly.com/javascript/error-bars/ -->
                    <div class="col-12" id="cases_survival_ai1" style="display: none;" >
                        <p> The plots below show the Kaplan meyer survival plot in each subgroup </p>
                        <div id="plot_cases_survival_ai1" class="mt-3 row justify-content-start"  >  </div>
                    </div>
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

function render_clinical_filter_projs_area_ai1(){
    project_filter_ai1.style.display = 'none';

    let domid_container = "filters_area_ai1";

    projects = [ "TCGA-ACC",  "TCGA-BLCA",  "TCGA-BRCA",  "TCGA-CESC",  "TCGA-CHOL",  "TCGA-COAD",  "TCGA-DLBC",  "TCGA-ESCA",  "TCGA-GBM",  "TCGA-HNSC",  "TCGA-KICH",  "TCGA-KIRC",  "TCGA-KIRP",  "TCGA-LAML",  "TCGA-LGG",  "TCGA-LIHC",  "TCGA-LUAD",  "TCGA-LUSC",  "TCGA-MESO",  "TCGA-OV",  "TCGA-PAAD",  "TCGA-PCPG",  "TCGA-PRAD",  "TCGA-READ",  "TCGA-SARC",  "TCGA-SKCM",  "TCGA-STAD",  "TCGA-TGCT",  "TCGA-THCA",  "TCGA-THYM",  "TCGA-UCEC",  "TCGA-UCS",  "TCGA-UVM" ]
    let domid_target = "project_ai1";
    let label = "Project";
    let options = projects;
    let selected = options[0];
    let onchange = ""
    fill_select( label, options, domid_target, domid_container, selected, onchange);

    domid_target = "dimension_ai1";
    label = "Demographic Variables";
    options = ["all", "race", "gender", "ethnicity"];
    selected = options[0];
    onchange = ""
    fill_select( label, options, domid_target, domid_container, selected, onchange);
}

function _render_pie_hist_stage(dat_cases){
    let sel_var = select_dimension_ai1.value;
    let container = "plot_cases_hist_stage_ai1";

    let targets = ['pathologic_stage'];
    let tmp = { 'all': {  } };
    for(let d of dat_cases){
        let v = d[sel_var];
        if( ! Object.keys(cnt).includes(v) ){
            tmp[v] = {};
        }
        for(let t of targets){
            t = d[t]

            if( ! Object.keys(tmp[v]).includes(t) ){
                tmp[v][t] = 0;
            }
            tmp[v][t] += 1

            if( ! Object.keys(tmp['all']).includes(t) ){
                tmp['all'][t] = 0;
            }
            tmp['all'][t] += 1
        }
    }

    let keys = Object.keys(tmp);
    let div_ids = [];
    let htmls = "";
    keys.forEach( (it) => {
        let _id = it.replaceAll(' ','_');

        div_ids.push(`sec1_plot_${_id}_ai1`);
        htmls += `
            <div  id = "pie_${_id}_ai1" class="col-4" >

            </div>
        `;
    });
    document.getElementById(container).innerHTML = htmls;
    
    keys.forEach( (it) => {
        let _id = it.replaceAll(' ','_');

        let itlay = pie_layout;
        itlay["title"] = { "text": `Subgroup ${ _id }` };
        let labels = Object.keys( tmp[it] );
        let values = Object.values( tmp[it] );
        let pldata = [{ "values": values, "labels": labels, "type": "pie" }];
        Plotly.newPlot( `pie_${_id}_ai1`, pldata, itlay, config);
    });

}

function _render_prescribed_drugs(dat_cases){
    let sel_var = select_dimension_ai1.value;
    let container = "plot_cases_prescribed_drugs_ai1";

    let targets = ['name'];
    let tmp = { 'all': { } };
    for(let d of dat_cases){
        let v = d[sel_var];
        if( ! Object.keys(cnt).includes(v) ){
            tmp[v] = {};
        }
        for(let drg of d['drug_details']){
            t = drg['name']
            if( ! Object.keys(tmp[v]).includes(t) ){
                tmp[v][t] = 0;
            }
            tmp[v][t] += 1

            if( ! Object.keys(tmp['all']).includes(t) ){
                tmp['all'][t] = 0;
            }
            tmp['all'][t] += 1
        }
    }
}

function _render_distribution_plots(dat_cases){
    let sel_var = select_dimension_ai1.value;
    let container = "plot_cases_dist_ai1";

    let all_targets = ['dosis_drugs', 'dosis_radiation', 'duration_treat_drug', 'duration_treat_radiation', 'age_at_initial_pathologic_diagnosis', 'qty_follow_ups', 'qty_drugs', 'qty_radiation'];
    let targets = ['age_at_initial_pathologic_diagnosis', 'qty_follow_ups', 'qty_drugs', 'qty_radiation'];
    let tmp = { 'all': { } };
    for(let d of dat_cases){
        let v = d[sel_var];
        if( ! Object.keys(cnt).includes(v) ){
            tmp[v] = { };
            for(let t of all_targets){
                tmp[v][t] = [];
                tmp['all'][t] = [];
            }
        }

        for(let t of targets){
            tmp[v][t].push( parseInt(d[sel_var][t]) );
            tmp['all'][t].push( parseInt(d[sel_var][t]) );
        }

        for(let drg of d['drug_details']){
            v1 = drg['duration'];
            v2 = drg['dosis'];
            if(v1 != null){
                tmp[v]['duration_treat_drug'].push( parseInt( v1 ) );
                tmp['all']['duration_treat_drug'].push( parseInt( v1 ) );
            }
            if(v2 != null){
                tmp[v]['dosis_drugs'].push( parseInt( v2 ) );
                tmp['all']['dosis_drugs'].push( parseInt( v2 ) );
            }
        }

        for(let rad of d['radiation_details']){
            v1 = rad['duration'];
            v2 = rad['dosis'];
            if(v1 != null){
                tmp[v]['duration_treat_radiation'].push( parseInt( v1 ) );
                tmp['all']['duration_treat_radiation'].push( parseInt( v1 ) );
            }
            if(v2 != null){
                tmp[v]['dosis_radiation'].push( parseInt( v2 ) );
                tmp['all']['dosis_radiation'].push( parseInt( v2 ) );
            }
        }
    }
}

function _render_km_survival_plots(dat_surv){
    let sel_var = select_dimension_ai1.value;
    let container = "plot_cases_survival_ai1";

}

function perform_render_stratification_analysis(){
    let selected_proj = select_project_ai1.value;

    // pie plot stage and horizontal bar of drugs per

    obj_cov.get_clinical_stratification_survival(selected_proj).then( (dat) => {
        
        tissue_ai1.innerHTML = dat.cases[0].tumor_tissue_site;
        hist_type_ai1.innerHTML = dat.cases[0].histological_type;

        let info = [];
        info.push(`<b>All samples</b>: ${dat.cases.length}`);
        let sel_var = select_dimension_ai1.value;
        if( sel_var != "all" ){
            let cnt = {};
            for(let d of dat.cases){
                let v = d[sel_var];
                if( ! Object.keys(cnt).includes(v) ){
                    cnt[v] = 0;
                }
                cnt[v] += 1
            }
            for(let v of Object.keys(cnt) ){
                info.push(`<b>${v}</b>: ${cnt[v]}`);
            }
        }
        samples_info_ai1.innerHTML = info.join(' | ')

        /*
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
        */
    } );
}

function setup_current_project_ai1(p){
    analysis_current_ai1.style.display = 'none';

    obj_ai1.current_project = obj_cov.projects.filter( x => x.project_id == p )[0];
    document.getElementById('proj_name_ai1').innerHTML = obj_cov.current_project.name;

    selected_proj.innerHTML = `Selected project: ${p}`;

    perform_render_stratification_analysis();

    analysis_current_ai1.style.display = '';  
}

let init_case_expai1 = async () => {
    document.getElementById('notice_ai1').innerHTML = 'Loading...';
      
    obj_ai1 = new GdcExplorerLib( );
    await obj_ai1.get_projects();
    await obj_ai1.get_all_project_summary( obj_ai1.pids );
    
    //render_filter_projs_area_ai1();
    render_clinical_filter_projs_area_ai1();

    document.getElementById('notice_ai1').innerHTML = '';
}

init_case_expai1().then( async v => { console.log("Initialized!"); } );


