class ExpAiDisparitySNVSummary extends HTMLElement {
  constructor() {
    super();
  }

  connectedCallback() {
    this.innerHTML = `
      <p id = "notice_ai2" > Loading... </p>

      <section id="expai2_analysis" class="mt-3">
            <p>
                Exploratory analysis concerning mutation annotations with group stratification, acccording to a project and a target demographic variable
            </p>
            
            <div id="project_filter_ai2" style="display: none;" > 
                <h4> Filter the project you want to explore </h4>

                <div class="row g-2 align-items-center" id="filters_area_ai2" >


                </div>

                <div class="mt-3 col-12" >
                    <button type="button" class="ml-3 btn btn-primary" id="go_filter" onclick="setup_current_project_ai2()" > Analyze </button>
                </div>
            </div>

            <hr />

            <div class="row justify-content-center mt-3"  >
                <div class="col-12" id="instructions_ai2" style="display: none;" >
                    <p>
                        Choose a project (clicking on the blue bar) to load details about the information available for the cases.
                    </p>
                </div>

                <div class="col-12" id="analysis_current_ai2" style="display: none;" >
                    <h4 class="mt-3" id="selected_proj_ai2" >  </h4>
                    <p>
                        <span style = "font-weight: bold;" > Name: </span> 
                        <span id = "proj_name_ai2" >  </span> 
                    </p>

                    <p>
                        <span style = "font-weight: bold;" > Tumor tissue site: </span> 
                        <span id = "tissue_ai2" >  </span> 
                    </p>

                    <p>
                        <span style = "font-weight: bold;" > Histological type: </span> 
                        <span id = "hist_type_ai2" >  </span> 
                    </p>

                    <p>
                        <span style = "font-weight: bold;" > Subgroups count: </span> 
                        <span id = "samples_info_ai2" >  </span> 
                    </p>
                    <p style="font-weight: bold;" > * the subgroups not represented in the next visualization sections did not have enough annotations in the cases clinical data</p>

                    <!-- plot mutation features grouped by subgroup -->
                    <div class="col-12" id="mut_features_ai2" style="display: none;" >
                        <h4> The plots below show the ranked frequency of mutation annotations found per subgroup</h4>
                        <div id="plot_mut_features_ai2" class="mt-3 row justify-content-start"  >  </div>
                    </div>

                    <!-- plot in_common and exclusive mutations per subgroup -->
                    <div class="col-12" id="intersection_exclusive_ai2" style="display: none;" >
                        <h4> The plots below show the ranked counts of mutations in common among the subgroups and the exclusive ones in each of these subgroups </h4>
                        <div id="plot_intersection_exclusive_ai2" class="mt-3 row justify-content-start"  >  </div>
                    </div>
        
                </div>

            </div>
        </section>
        
    `;
  }
}
customElements.define('expai2-component', ExpAiDisparitySNVSummary);

/*
Page organization:
sec 1 -  most frequent mutations , intersection and exclusive subgroups ranking - same for gene  - if all is not selected hide intersection and diff
sec 2 - distribution of impact, consequence, pathogenicity , grouped by subgroup
sec 3 - demographic arrangement of cases. given a project (select) there is also a select to filter the view by race, gener, ethnicity, age range (baby, child, teen, adult, elderly)
*/

let obj_ai2 = {};

let config_ai2 = { "responsive": true };
let pie_layout_ai2 = { "width": 500, "height": 500 };

// Section 0 - projects filtering
// obj_ai2.projects.filter( e => e.disease_type.map( e => e.toLowerCase() ).filter( it => it.indexOf('epith')!=-1 ) )

function _fill_data_category_ai2(){
    let p = select_project_ai2.value;
    let allowed = new Set( Object.keys(obj_ai2.formats_datCategory) )

    let domid_container = "filters_area_ai2";

    let domid_target = "data_category_ai2";
    let label = "Data Category";
    let options = new Set( obj_ai2.projects_summary[p].data_categories.map( e => e.data_category.toLowerCase() ) );
    options = options.intersection(allowed); // Filter data categories that have open parsable (not platform-specific raw output as idat) data
    options = Array.from( options );
    let selected = options[0];
    fill_select( label, options, domid_target, domid_container, selected, null);
}

function render_filter_projs_area_ai2(){
    project_filter_ai2.style.display = 'none';

    let domid_container = "filters_area_ai2";

    let domid_target = "project_ai2";
    let label = "Project";
    let options = Object.keys(obj_ai2.projects_summary);
    
    let selected = options[0];
    let onchange = "_fill_data_category()"
    fill_select( label, options, domid_target, domid_container, selected, onchange);

    _fill_data_category_ai2();

    project_filter_ai2.style.display = '';
}

function render_clinical_filter_projs_area_ai2(){
    project_filter_ai2.style.display = 'none';

    let domid_container = "filters_area_ai2";

    // Filtered by those that have general_snprs_info.json
    projects = [ "TCGA-ACC", "TCGA-BLCA", "TCGA-BRCA", "TCGA-CESC", "TCGA-CHOL", "TCGA-COAD", "TCGA-DLBC", "TCGA-ESCA", "TCGA-GBM", "TCGA-HNSC", "TCGA-KICH", "TCGA-KIRC", "TCGA-KIRP", "TCGA-LAML", "TCGA-LGG", "TCGA-LIHC", "TCGA-LUAD", "TCGA-LUSC", "TCGA-MESO", "TCGA-OV", "TCGA-PAAD", "TCGA-PCPG", "TCGA-PRAD", "TCGA-READ", "TCGA-SARC", "TCGA-SKCM", "TCGA-STAD", "TCGA-TGCT", "TCGA-THCA", "TCGA-THYM", "TCGA-UCEC", "TCGA-UCS", "TCGA-UVM" ]
    let domid_target = "project_ai2";
    let label = "Project";
    let options = projects;
    let selected = options[0];
    let onchange = ""
    fill_select( label, options, domid_target, domid_container, selected, onchange);

    domid_target = "dimension_ai2";
    label = "Demographic Variables";
    options = ["all", "race", "gender", "ethnicity"];
    selected = options[0];
    onchange = ""
    fill_select( label, options, domid_target, domid_container, selected, onchange);

    project_filter_ai2.style.display = '';
}

function _render_features_mutations( dat ){
    let container = "plot_mut_features_ai2";
    let sel_var = select_dimension_ai2.value;
    dat = dat[`by_${sel_var}`];
    subgroups = Object.keys(dat);

    let mp_name_info = { "impact": "Impact", "consequence": "Consequence", "clin_sig": "Pathogenicity", "hugo_symbol": "Gene" };
    let plot_types = Object.keys(mp_name_info);

    let layout = { title: { text: "Mutation features by " }, xaxis: { title: { text: 'Feature' } }, yaxis: { title: { text: 'Values' }, zeroline: false }, barmode: 'group' };

    plot_types.forEach( (pt) => {
        let ptname = pt;
        if( Object.keys(mp_name_info).includes(pt) ){
            ptname = mp_name_info[pt];
        }
        layout["title"] = { text: `Mutation features by ${ptname}` };

        let htmls = "";
        htmls += `
            <div  id = "bar_mut_feature_by_${pt}_ai1" class="col-6" >

            </div>
        `;
        document.getElementById(container).innerHTML += htmls;
        let pldata = [];
        subgroups.forEach( (it) => {
            let tmp = dat[it][pt];
            tmp = _sort_dict_by_values(tmp);

            let x = Object.keys(tmp);
            let y = Object.values(tmp);
            let ridx = x.indexOf('nan');
            if(ridx != -1){
                x.splice(ridx, 1);
                y.splice(ridx, 1);
            }
            if( x.length > 20 ){
                x = x.slice(0,20);
                y = y.slice(0,20);
            }
            
            let obj = { type: 'bar', x: x, y: y, name: it };
            pldata.push(obj)
        });

        Plotly.newPlot( `bar_mut_feature_by_${pt}_ai1`, pldata, layout, config);
    })
    

    document.getElementById("mut_features_ai2").style.display = '';
}

function perform_render_stratification_analysis_ai2(){
    let selected_proj = select_project_ai2.value;

    // pie plot stage and horizontal bar of drugs per

    obj_ai2.load_cases_mutation_annotations(selected_proj).then( (dat) => {

        tissue_ai2.innerHTML = dat.cases[0].tumor_tissue_site;
        hist_type_ai2.innerHTML = dat.cases[0].histological_type;

        let info = [];
        info.push(`<b>All samples</b>: ${dat.cases.length}`);
        let sel_var = select_dimension_ai2.value;
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
        samples_info_ai2.innerHTML = info.join(' | ')

        _render_features_mutations(dat.ann_mutations);

        /*
         _render_pie_hist_stage(dat.cases);
         _render_prescribed_drugs(dat.cases);
         _render_distribution_plots(dat.cases);
         _render_km_survival_plots(dat.survival);
         */
    } );
}

function setup_current_project_ai2(){
    analysis_current_ai2.style.display = 'none';
    let p = select_project_ai2.value;

    obj_ai2.current_project = obj_ai2.projects.filter( x => x.project_id == p )[0];
    document.getElementById('proj_name_ai2').innerHTML = obj_ai2.current_project.name;

    selected_proj.innerHTML = `Selected project: ${p}`;

    perform_render_stratification_analysis_ai2();

    analysis_current_ai2.style.display = '';  
}

let init_case_expai2 = async () => {
    document.getElementById('notice_ai2').innerHTML = 'Loading...';
      
    obj_ai2 = new GdcExplorerLib( );
    await obj_ai2.get_projects();
    await obj_ai2.get_all_project_summary( obj_ai2.pids );
    
    render_clinical_filter_projs_area_ai2();

    document.getElementById('notice_ai2').innerHTML = '';
}

init_case_expai2().then( async v => { console.log("Initialized!"); } );


