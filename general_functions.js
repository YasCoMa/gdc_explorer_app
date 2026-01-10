function _capitalize(vs){
    let nv = vs[0].toUpperCase() + vs.slice(1);
    nv = nv.replaceAll('_',' ');

    return nv;
}

function _get_exclusive_values(target_key, obj, type_value = "dict"){
	let tvalues = obj[target_key];
	if( type_value == "dict" ){
		tvalues = Object.keys(obj[target_key]);
	}
	let to_remove = new Set();

	let values = [];
	for(let k of Object.keys(obj) ){
		if( k != target_key ){
			values = obj[k];
			if( type_value == "dict" ){
				values = Object.keys(obj[k]);
			}
			for( let v of values){
				if( tvalues.includes(v) ){
					to_remove.add(v);
				}
			}
		}
	}
	tvalues = new Set(tvalues);
	tvalues = tvalues.difference(to_remove);

	let exclusive = Array.from(tvalues);

	return exclusive;
}

function _get_intersection_values(obj, type_value = "dict"){
	let keys = Object.keys(obj);
	let tvalues = obj[ keys[0] ];
	if( type_value == "dict" ){
		tvalues = Object.keys(obj[ keys[0] ]);
	}
	let inter = new Set(tvalues);

	let values = [];
	for(let k of Object.keys(obj) ){
		if( k != target_key ){
			values = obj[k];
			if( type_value == "dict" ){
				values = Object.keys(obj[k]);
			}
			let nv = new Set(values);
			inter = inter.intersection(nv);
		}
	}

	let in_common = Array.from(inter);

	return in_common;
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
