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
                 The approved evidence described in the release launched in June 2026 was filtered for that which is “Predictive” and that has the “Supports” direction.
            </p>
            <button type="button" class="mt-3 btn btn-secondary" onclick="if( read_more_ai3.style.display == 'none' ){ read_more_ai3.style.display = ''; } else { read_more_ai3.style.display = 'none' }" > Read more </button>

            <div id = "read_more_ai3" style = "margin-top: 10px; display: none;" >
                <h5>Context</h5>
                <p>
                The nucleus of our cells contains the genome, a series of letters (nucleotides) that encodes each individual unique genetic code with instructions to produce the proteins and how each body tissue and organ must behave. Sometimes, there are a pile of letter changes (mutations) in specific regions of the genome that are not properly repaired by the DNA proofreading mechanism, or they are severely accumulated, and this may lead to the production of proteins that do not work correctly and later induce the uncontrolled growth of undesired cells.
            <br />
                In the first stages of cancer (I, II), blood tests of predictive and diagnostic genes may detect it very early, or by image and biopsy exams (to make sure that the tissue and cell have cancerigenous properties). When the mass is small, the prognosis and treatment are less aggressive, and surgery and radiotherapy are able to solve it. However, when it is not detected at this stage, the cells continue to divide and reproduce, and at a certain point, the local resources are not able to nourish them, so they need to create their own mechanism to receive the nutrients, so they create their own circulatory system with new blood vessels (angiogenesis).
            <br />
                In stage III, the tumor cells start to migrate and spread to local boundary tissues and local lymph nodes. In stage IV, tumor circulating cells are generated, capable of entering the bloodstream and lymphatic system in a process named metastasis. These cells may migrate to any other part of the body and create new colonies, with new characteristics that may turn the disease even more complex to treat depending on the protein functional profile. Up to stage III, it is still possible to combine surgical intervention with radio and chemo; however, in the metastatic stage, conventional treatments are hard to achieve complete remission.
            <br />
                New lines of treatment are being derived through precision oncology, where the clinical decisions of treatments (either drugs, vaccines, or surgical interventions) are guided according to the whole set of genetic and protein alterations that the patients have. There are specific hotspots in known predictive and biomarker genes for which some drugs were approved, targeting the proteins with the main driver structural phenotype and causing their inhibition, which is able to affect downstream signaling cascades and lead to tumor reduction and remission. 
            <br />
                Apart from the mutations in genes that eventually become potential targets for drug binding, there are other sets of mutations in genes that functionally act as enzymes and are responsible for degrading chemical compounds (drugs) and releasing the desired substrates in the bloodstream. The drugs that we take follow a standard path (ADMET) that begins with its Absortion (the way it is administered / enters the organism), following with their distribution (how the substance will traverse the bloodstream or cross the blood-brain barrier and reach the location that it has to produce the effect, binding to a target); once in the correct location and reaching the final structural conformation, it has to be Metabolized (broken in pieces) by the its specific set of enzymes to release the substrates that are needed to activate the planned regulatory pathways; the remaining substrates of the metabolism phase must be Eliminated from the body, usually being filtered in the kidneys and terminating as part of the urine components; lastly, but extremely important, it is the evaluation of the drug's toxicity, assessing the possible adverse reactions, adjusting of the dosage according to the patient response, how dangerous the expected effects of its metabolism products are, and which organs and blood markers should be monitored to track and account for severe systemic harm. 
            <br />
                The mutations in these last genes may affect the speed of a certain drug degradation, and depending on the type of drug and the dosage, it may accumulate and become toxic for people carrying the mutation. Precision oncology is also important in this aspect, to provide prescription guidance, adjust the dosage of first-line treatment drugs for a condition, or even recommend other treatment options when the deficiency is even more severe.

            <br />
            
            <b>Sources</b>: 
            <ul>
                <li>https://cancer.ca/en/cancer-information/what-is-cancer/how-cancer-starts-grows-and-spreads</li>
                <li>https://www.sciencedirect.com/topics/pharmacology-toxicology-and-pharmaceutical-science/absorption-distribution-metabolism-excretion-toxicity</li>
            </ul>

                </p>
            </div>
            
            <div id="project_filter_ai3" style="margin-top: 15px; display: none;" > 
                <h4> Filter the project you want to explore </h4>

                <div class="row g-2 align-items-center" id="filters_area_ai3" >


                </div>

                <div class="mt-3 col-12" >
                    <button type="button" class="ml-3 btn btn-primary" id="go_filter" onclick="setup_current_project_ai3()" > Analyze </button>
                </div>
            </div>

            <hr />

            <div class="row justify-content-center mt-3"  >
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

                    <div class="col-12" id="analysis_clinpgx_ann_ai3" style="display: none;" >
                        <h4 class="mt-3" > Annotations from ClinPGx database for the cases </h4>

                        <div class="col-12" id="cases_summary_table_ai3_clinpgx" style = "display: none;"  >
                            <h5>  <b> Summary annotations </b> </h5>

                            <div class="scroll-container" >
                                <table class="table table-striped scroll-container" id="clinpgx_summary_datatab" >
                                
                                </table>
                            </div>
                        </div>

                        <div class="col-12 mt-3" id="cases_variant_table_ai3_clinpgx" style = "display: none;"  >
                            <h5>  <b> Variant annotations </b> </h5>

                            <div class="scroll-container" >
                                <table class="table table-striped" id="clinpgx_variant_datatab" >
                                    
                                </table>
                            </div>
                        </div>

                        <hr />
                    </div>

                    <div class="col-12 mt-2"  >
                        <h4 class="mt-3" > Explore mutation information from CIVICdb </h4>
                        
                        <p id="instructions_ai3">
                            Choose a mutation (clicking on the blue bar) to load details about the information available for the cases.
                        </p>

                        <div  class="mt-3 row justify-content-center"  >
                            <div id="plot_mutation_case_count" class="col-12 " >  </div>
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
        
        <style>
        .scroll-container {
            max-height: 500px;        /* 1. You must define a fixed limit */
            overflow-y: auto;     /* 2. Adds a vertical scrollbar ONLY when content overflows */
            overflow-x: hidden;   /* 3. Optional: Prevents accidental horizontal scrolling */
        }
        </style>
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
    let onchange = "setup_current_project_ai3()"
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

function _render_cases_table_for_summary_ann(){
    document.getElementById("clinpgx_summary_datatab").innerHTML = "";
    document.getElementById("cases_summary_table_ai3_clinpgx").style.display = 'none';

    let htmls = `
        <tr>
            <th> Case </th>
            <th> Allele </th>
            <th> Function </th>
            <th> Diseases </th>
            <th> Drugs </th>
            <th> Annotation </th>
        </tr>
    `;

    let dat = obj_ai3.annotations_pgx;
    let cases = Object.keys(dat).filter( k => Object.keys(dat[k]['summary_clingx']).length > 0 );
    if( cases.length > 0 ){
        for( let c of cases ){
            let d = dat[c]["summary_clingx"];
            let alleles = Object.keys(d);
            for(let allele of alleles){
                let function_al = d[allele]['function'];

                let annots = d[allele]['anns'][0];
                for(let ann of annots){
                    let sentence = ann["ann"];
                    sentence = `<p style="font-size: 0.8rem !important" > ${sentence} </p>`;
                    
                    let diseases = ann['diseases'].map( k => `<span class="badge bg-primary text-light"> ${k} </span>` ).join(', ');

                    let drugs = ann['drugs'].map( k => `<span class="badge  bg-info text-dark"> ${k} </span>` ).join(', ');

                    htmls += `
                        <tr>
                            <td> ${ c } </td>
                            <td> ${ allele } </td>
                            <td> ${ function_al } </td>
                            <td> ${ diseases } </td>
                            <td> ${ drugs } </td>
                            <td> ${ sentence } </td>
                        </tr>
                    `;
                }

            }
        }
        document.getElementById("clinpgx_summary_datatab").innerHTML = htmls;

        document.getElementById("cases_summary_table_ai3_clinpgx").style.display = '';
    }
}

function _render_cases_table_for_variant_ann(){
    document.getElementById("clinpgx_variant_datatab").innerHTML = "";
    document.getElementById("cases_variant_table_ai3_clinpgx").style.display = 'none';

    let htmls = `
        <tr>
            <th> Case </th>
            <th> Variant </th>
            <th> Genotype </th>
            <th> Phenotype Category </th>
            <th> Drugs </th>
            <th> Annotation </th>
        </tr>
    `;

    let dat = obj_ai3.annotations_pgx;
    let cases = Object.keys(dat).filter( k => dat[k]['variant_clingx'].length > 0 );
    if( cases.length > 0 ){
        for( let c of cases ){
            let annots = dat[c]["variant_clingx"];
            for(let ann of annots){
                let gene = ann["gene"];
                let rsid = ann["rsid"];
                let variant = `${rsid} in gene ${gene}`;

                let genotype = ann["genotype"][0];

                let sentence = ann["ann"];
                sentence = `<p style="font-size: 0.8rem !important" > ${sentence} </p>`;
                
                let phenoCats = ann['phenoCats'].map( k => `<span class="badge bg-primary text-light"> ${k} </span>` ).join(', ');

                let drugs = ann['drugs'].map( k => `<span class="badge bg-info text-dark"> ${k} </span>` ).join(', ');

                htmls += `
                    <tr>
                        <td> ${ c } </td>
                        <td> ${ variant } </td>
                        <td> ${ genotype } </td>
                        <td> ${ phenoCats } </td>
                        <td> ${ drugs } </td>
                        <td> ${ sentence } </td>
                    </tr>
                `;
            }
        }
        document.getElementById("clinpgx_variant_datatab").innerHTML = htmls;

        document.getElementById("cases_variant_table_ai3_clinpgx").style.display = '';
    }
}

function render_clinpgx_annotations_details(){
    _render_cases_table_for_summary_ann();
    _render_cases_table_for_variant_ann();

    analysis_clinpgx_ann_ai3.style.display = '';
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

    obj_ai3.load_pgx_mutation_data(selected_proj).then( (resp) => {
        let dat = resp.result;

        tissue_ai3.innerHTML = dat.tissue_site;
        hist_type_ai3.innerHTML = dat.hist_type;
        obj_ai3.data_pgx = dat;
        obj_ai3.annotations_pgx = resp.annotations;

        _render_overview_count_mutations(dat.mutation_count);
        let all_muts = Object.keys(dat.mutation_count);
        render_mutation_details(all_muts[0]);

        clinpgx_keys = Object.keys(obj_ai3.annotations_pgx);
        if(clinpgx_keys.length > 0){
            render_clinpgx_annotations_details();
        }

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


