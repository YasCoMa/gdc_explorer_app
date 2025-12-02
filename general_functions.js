function _capitalize(vs){
    let nv = vs[0].toUpperCase() + vs.slice(1);
    nv = nv.replaceAll('_',' ');

    return nv;
}

const intersect = (setA, setB, ...args) => {
  const result = new Set([...setA].filter((i) => setB.has(i)))
  if (args.length === 0) return result
  return intersect(result, args.shift(), ...args)
}

function _sum(arr){
    return arr.reduce(function (a, b) { return a + b; }, 0);
}

function _make_assoc_array(e, header){
	let d = {};
	let i = 0;
	for ( let h of header){
		d[h] = e[i];
		i += 1;
	}
	
	return d;
}

function fill_select( label, options, domid_target, domid_container, selected, function_change){
    let fadd = '';
    if( function_change ){
        fadd = ` onChange="${function_change}" `
    }
    
	let ops = "";
	options.forEach( e => {
		let add = "";
		if( e == selected ){
			add += "selected";
		}
		ops += `<option value = "${ e }" ${ add } > ${ e } </option>` 
	} );

	let eldom = document.getElementById(`select_${domid_target}`);
	if( ! eldom ){
		let htmls = `
			<div class="col-auto">
		        <label class="form-label" > ${ label }:</label>
				<select class="form-select" id="select_${domid_target}" ${fadd} >
					${ ops }
				</select>
			</div>
		`;
		document.getElementById(domid_container).innerHTML += htmls;
	}
	else{
		eldom.innerHTML = ops;
	}
}

function generate_tags(options, domid_container){
	let htmls = "";
	options.forEach( e => {
		htmls += `<span class="badge rounded-pill bg-primary mx-2"> ${ e } </span>` 
	} );

	document.getElementById(domid_container).innerHTML = htmls;
}
