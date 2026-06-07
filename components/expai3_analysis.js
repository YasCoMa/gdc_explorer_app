class ExpAiCivicSNVSummary extends HTMLElement {
  constructor() {
    super();
  }

  connectedCallback() {
    this.innerHTML = `
      <p id = "notice_ai3" > Loading... </p>

      <section id="expai3_analysis" class="mt-3">
            <p>
                Exploratory analysis concerning mutations found in the TCGA project cases that have some evidence of association with drugs in <a href="https://civicdb.org/releases/main" target="_blank" > CIVIC database</a>.
                <br />
                The approved evidences described in the release launched on June, 2026 were filtered for those that are Predictive and that has the Supports direction.
             </p>
            
            <div id="project_filter_ai3" style="display: none;" > 
                <h4> Filter the project you want to explore </h4>

                <div class="row g-2 align-items-center" id="filters_area_ai3" >


                </div>

                <div class="mt-3 col-12" >
                    <button type="button" class="ml-3 btn btn-primary" id="go_filter" onclick="setup_current_project_ai3()" > Analyze </button>
                </div>
            </div>

            <hr />

            <div class="row justify-content-center mt-3"  >
                <div class="col-12" id="instructions_ai3" style="display: none;" >
                    <p>
                        Choose a mutation (clicking on the blue bar) to load details about the information available for the cases.
                    </p>
                </div>

                <div class="col-12" id="analysis_current_ai3" style="display: none;"  >
                    <h4 class="mt-3" id="selected_proj_ai3" >  </h4>
                    <p>
                        <span style = "font-weight: bold;" > Name: </span> 
                        <span id = "proj_name_ai3" >  </span> 
                    </p>

                    <p>
                        <span style = "font-weight: bold;" > Tumor tissue site: </span> 
                        <span id = "tissue_ai3" >  </span> 
                    </p>

                    <p>
                        <span style = "font-weight: bold;" > Histological type: </span> 
                        <span id = "hist_type_ai3" >  </span> 
                    </p>
                    
                    <div class="col-12"  >
                        <p>
                            Choose a mutation (clicking on the blue bar) to load details about the information available for the cases.
                        </p>

                        <div  class="mt-3 row justify-content-center"  >
                            <div id="plot_mutation_case_count" class="col-12 " >  </div>
                        </div>
                    </div>

                    <div class="col-12" id="analysis_current_mutation_ai3" style="display: none;" >
                        <h4 class="mt-3" id="selected_mut_ai3" >  </h4>

                        <div class="col-12" id="mut_description_ai3" style="display: none;" >
                        </div>

                        <!-- plot demographic information of cases -->
                        <div class="col-12" id="cases_overview_plot_ai3" style = "display: none;" >
                            <p> The plots below show the distribution of cases across demographic information. Some cases do not have some annotations such as the disease stage. </p>
                            <div id="plot_cases_overview_plot_ai3" class="mt-3 row justify-content-start"  >  </div>
                        </div>

                        <!-- table cases details -->
                        <div class="col-12" id="cases_table_ai3" style = "display: none;"  >
                            <p> The table below show the information of cases that have this mutation. The information with * represent items that were matched in CIVIC database. </p>
                            
                            <div  class="mt-3 row justify-content-start"  >
                                <table class="table table-striped" id="datatab" >
                                
                                </table>
                            </div>
                        </div>
                    </div>
        
                </div>

            </div>
        </section>

        <div id="full-mutation-list" class="modal" tabindex = "-1" >
            <div class="modal-dialog modal-dialog-centered modal-dialog-scrollable" >
                <div class="modal-content">

                  <div class="modal-header">
                    <h1 class="modal-title fs-5" id="staticBackdropLabel">
                        Full list of mutations
                    </h1>

                    <button type="button" class="btn-close" data-bs-dismiss="modal" aria-label="Close"></button>
                  </div>

                  <div class="modal-body" id="full_muts" >
                        
                  </div>

                  <div class="modal-footer">
                    <button type="button" class="btn btn-secondary" data-bs-dismiss="modal">Close</button>
                  </div>

                </div>
            </div>
        </div>
        
    `;
  }
}
customElements.define('expai3-component', ExpAiCivicSNVSummary);

/*
Page organization:
sec 1 -  most frequent mutations , intersection and exclusive subgroups ranking - same for gene  - if all is not selected hide intersection and diff
sec 2 - distribution of impact, consequence, pathogenicity , grouped by subgroup
sec 3 - demographic arrangement of cases. given a project (select) there is also a select to filter the view by race, gener, ethnicity, age range (baby, child, teen, adult, elderly)
*/

let obj_ai3 = {};

let config_ai3 = { "responsive": true };
let pie_layout_ai3 = { "width": 500, "height": 500 };

// Section 0 - projects filtering
// obj_ai3.projects.filter( e => e.disease_type.map( e => e.toLowerCase() ).filter( it => it.indexOf('epith')!=-1 ) )

function render_clinical_filter_projs_area_ai3(){
    project_filter_ai3.style.display = 'none';

    let domid_container = "filters_area_ai3";

    // Filtered by those that have general_snprs_info.json
    projects = [ "TCGA-BLCA", "TCGA-BRCA", "TCGA-CESC", "TCGA-CHOL", "TCGA-COAD", "TCGA-DLBC", "TCGA-ESCA", "TCGA-GBM", "TCGA-HNSC", "TCGA-KIRC", "TCGA-LGG", "TCGA-LIHC", "TCGA-LUAD", "TCGA-LUSC", "TCGA-MESO", "TCGA-OV", "TCGA-PAAD", "TCGA-PRAD", "TCGA-READ", "TCGA-SARC", "TCGA-SKCM", "TCGA-STAD", "TCGA-TGCT", "TCGA-THCA", "TCGA-UCEC", "TCGA-UCS", "TCGA-UVM" ]
    let domid_target = "project_ai3";
    let label = "Project";
    let options = projects;
    let selected = options[0];
    let onchange = ""
    fill_select( label, options, domid_target, domid_container, selected, onchange);

    project_filter_ai3.style.display = '';
}

function _render_mutation_description(mut){
    let container = "mut_description_ai3";
    let htmls = "";
    let dat = obj_ai3.data_pgx.mutation_details[mut];
    let names = ["Chromosome", "Gene", "Variant consequence", "PolyPhen annotation", "Clinical meaning", "Associated diseases in CIVIC", "Associations to drugs"];
    let descriptors = ["chromosome", "gene", "consequence", "polyphen", "clinsig", "civic_diseases", "drug-variant_association"];
    let i = 0;
    for(let d of descriptors){
        let v = dat[d];
        v = v.replaceAll("Resistance", `<span style="color: #D92344;" > Resistance </span>` ).replaceAll("Sensitivity/Response", `<span style="color: #47AD40;" > Sensitivity/Response </span>` );
        let title = names[i];
        if( v!="nan" && v!="" && v!=null && v!=undefined ){
            htmls += `<p> <span style="font-weight: bold;" > ${title} </span> ${v} </p>`;
        }
        i+=1;
    }
    document.getElementById(container).innerHTML = htmls;
    mut_description_ai3.style.display = '';
}

function _render_pie_mut_overview(mut){
    let container = "plot_cases_overview_plot_ai3";

    let tmp = obj_ai3.data_pgx.mutation_count[mut];
    let keys = Object.keys(tmp);
    keys.splice( keys.indexOf("cases"), 1 )
    let htmls = "";
    keys.forEach( (it) => {
        let _id = it.replaceAll(' ','_');
        htmls += `
            <div  id = "pie_${_id}_ai3" class="col-4" >

            </div>
        `;
    });
    document.getElementById(container).innerHTML = htmls;
    
    let uniques_null = 0;
    keys.forEach( (it) => {
        let _id = it.replaceAll(' ','_');

        let itlay = pie_layout;
        itlay["title"] = { "text": `Descriptor ${ _id }` };
        let labels = Object.keys( tmp[it] );
        if( labels.length == 1 && labels[0]=='null' ){
            uniques_null += 1;
        }
        else{
            let values = Object.values( tmp[it] );
            let pldata = [{ "values": values, "labels": labels, "type": "pie" }];
            Plotly.newPlot( `pie_${_id}_ai3`, pldata, itlay, config);
        }
        console.log(it, labels)
    });
    if(uniques_null == keys.length){
        document.getElementById(container).innerHTML = "<p> There is no information available about the cases nvolving this mutation. </p>";
    }

    document.getElementById("cases_overview_plot_ai3").style.display = '';
}

function _fill_mutation_modal(case_id){
    let mutations = obj_ai3.data_pgx.cases_data[case_id]["allmuts"];

    let mlist = mutations.map( m => `<li> <span style="font-weight: bold;"> ${ m.split("-")[0] } </span> - <span style="font-weight: bold;"> ${ m.split("-")[1] } </span> </li>` );
    mlist = mlist.join('');
    let htmls = `
        <ul>
            ${ mlist }
        </ul>
    `;
    document.getElementById("full_muts").innerHTML = htmls;
    let options = { focus: true };
    let modal = new bootstrap.Modal('#full-mutation-list', options);
    modal.show();

}

function _get_interplay_drugs_response_for_case(pdrugs, muts){
    let tpdrugs = pdrugs.map( e => e.replace("*", "").toLowerCase() );

    let nd = [];
    let gone = [];
    for( let m of muts ){
        let tmp = [];
        let dr = obj_ai3.data_pgx.mutation_details[m]["drug-variant_association"].split(", ").map( e => e.split(" to ") );
        for( let d of dr ){
            if( tpdrugs.includes( d[1].toLowerCase() ) ){
                let color = "#47AD40";
                if(d[0].includes("Resist")){
                    color = "#D92344";
                }
                let ann = `${d[1]} (<span style="color: ${color};" > ${d[0]} </span> )`;
                if( ! nd.includes(ann) ){
                    nd.push( ann )
                    gone.push(d[1])
                    gone.push( d[1].toLowerCase() )
                }
            }
        }
    }

    for(let d of pdrugs){
        if( ! gone.includes(d.replace("*","")) && d!="nan" ){
            nd.push(d)
        }
    }

    nd = nd.join(", ");
    return nd;
}

function _render_cases_table_for_current_mut(mut){
    let htmls = `
        <tr>
            <th> Case </th>
            <th> Demographic information </th>
            <th> Chromosomes </th>
            <th> Genes </th>
            <th> Mutations </th>
            <th> Prescribed drugs </th>
        </tr>
    `;

    let dat = obj_ai3.data_pgx.cases_data;
    let cases = obj_ai3.data_pgx.mutation_details[mut]["cases"];
    for( let c of cases ){
        let d = dat[c];
        let democols = ["race", "gender", "ethnicity"];
        let demos = [];
        for( let col of democols){
            if( d[col] && d[col]!="nan" ){
                demos.push( `<span style="font-weight: bold;"> ${ _capitalize(col) } </span> ${d[col]}` );
            }
        }
        demos = demos.join(" | ");
        let mutations = d["allmuts"].join("|");

        let pdrugs = d["prescribed_drugs"];
        let muts = d["muts"];
        console.log(muts)
        let ann_drugs = _get_interplay_drugs_response_for_case(pdrugs, muts);

        htmls += `
        <tr>
            <td> ${ c } </td>
            <td> ${ demos } </td>
            <td> ${ d["chromosomes"] } </td>
            <td> ${ d["genes"] } </td>
            <td> 
                <span> [${ d["muts"].join(', ') }] (CIVIC) / ${ d["allmuts"].length } (all) </span>
                <a href="javascript:void(0)" onclick="_fill_mutation_modal('${c}');" > See all mutations </a>
            </td>
            <td> ${ ann_drugs } </td>
        </tr>
        `;

    }
    document.getElementById("datatab").innerHTML = htmls;

    document.getElementById("cases_table_ai3").style.display = '';
}

function render_mutation_details(mut){
    obj_ai3.current_mutation = mut;

    selected_mut_ai3.innerHTML = mut;
    _render_mutation_description(mut);
    _render_pie_mut_overview(mut);
    _render_cases_table_for_current_mut(mut);

    analysis_current_mutation_ai3.style.display = '';
    console.log(mut)
}

function _render_overview_count_mutations( dat ){
    let layout = { barmode: 'group', title: { text: "Cases count per Mutation" }, xaxis: { title: { text: 'Mutations' } }, yaxis: { title: { text: 'Cases Count' } } }
    let pl = [];
    let x = Object.keys(dat);
    let y = x.map( e => dat[e].cases );
    pl.push( { "name": "Cases", "x": x, "y": y, "type": "bar" } );
    
    let sec1plot = document.getElementById('plot_mutation_case_count');
    Plotly.newPlot( "plot_mutation_case_count", pl, layout, config_ai3 );

    sec1plot.on('plotly_click', function(data){
        let p = data.points[0].x;

        render_mutation_details(p);
        
    });
    
    document.getElementById("instructions_ai3").style.display = '';
}

function perform_render_stratification_analysis_ai3(){
    let selected_proj = select_project_ai3.value;

    // pie plot stage and horizontal bar of drugs per

    obj_ai3.load_pgx_mutation_data(selected_proj).then( (dat) => {
        tissue_ai3.innerHTML = dat.tissue_site;
        hist_type_ai3.innerHTML = dat.hist_type;
        obj_ai3.data_pgx = dat;

        _render_overview_count_mutations(dat.mutation_count);
        let all_muts = Object.keys(dat.mutation_count);
        render_mutation_details(all_muts[0]);

        /*
         _render_pie_hist_stage(dat.cases);
         _render_prescribed_drugs(dat.cases);
         _render_distribution_plots(dat.cases);
         _render_km_survival_plots(dat.survival);
         */
    } );
}

function setup_current_project_ai3(){
    analysis_current_ai3.style.display = 'none';
    let p = select_project_ai3.value;

    obj_ai3.current_project = obj_ai3.projects.filter( x => x.project_id == p )[0];
    document.getElementById('proj_name_ai3').innerHTML = obj_ai3.current_project.name;

    selected_proj_ai3.innerHTML = `Selected project: ${p}`;

    perform_render_stratification_analysis_ai3();

    analysis_current_ai3.style.display = '';  
}

let init_case_expai3 = async () => {
    document.getElementById('notice_ai3').innerHTML = 'Loading...';
      
    obj_ai3 = new GdcExplorerLib( );
    await obj_ai3.get_projects();
    await obj_ai3.get_all_project_summary( obj_ai3.pids );
    
    render_clinical_filter_projs_area_ai3();
    setup_current_project_ai3();

    document.getElementById('notice_ai3').innerHTML = '';
}

init_case_expai3().then( async v => { console.log("Initialized!"); } );


