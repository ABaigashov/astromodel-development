

// ------------------------------------------ //


const defaultTextInputCases = {
	class: 'form-control',
	autocomplete: 'off',
	type: 'text'
}

const defaultNumereticInputCases = {
	class: 'form-control',
	autocomplete: 'off',
	type: 'number',
	step: 'any'
}

// ------------------------------------------ //

function show404() {
	document.body.innerHTML = '<h1>Error 404: Page not found</h1>'
}

function localhost(htmlForm) {
	let fileToLoad = Array.from(document.getElementById('filedrop').files).filter((file) => (file.name == 'config.json'))[0]
	let fileReader = new FileReader()
	let problemName = fileToLoad.webkitRelativePath.split('/')[0]
	fileReader.onload = (event) => {
		let config = JSON.parse(event.target.result)
		document.body.removeChild(document.body.lastElementChild)
		document.getElementById('problem-form').removeAttribute('style')
		generateAllProblemForm(problemName, config)
	}
	fileReader.readAsText(fileToLoad, 'UTF-8')
}

function callGenerate() {
	if (window.location.protocol === 'file:') {
		let div = document.createElement('div')
		document.body.appendChild(div)
		document.getElementById('problem-form').setAttribute('style', 'display:none;')
		div.innerHTML = '<input type="file" onchange="localhost()" style="display:none;" id="filedrop" webkitdirectory multiple /><label for="filedrop" style="position:absolute;width:200px;height:100px;margin:auto;left:0;top:0;right:0;bottom:0;display:flex;justify-content:center;align-items:center;background-color:#bbb;border-style:dashed">Drop problem directory</label>'
	} else {
		let problemName = window.location.search.substring(1)
		$.ajax({
			url: window.location.protocol + '//' + window.location.host + '/construct/configs/' + problemName + '.json',
			type: 'GET',
			dataType: 'json',
			contentType: 'application/json; charset=utf-8',
			success: function(config) {
				generateAllProblemForm(problemName, config)
			},
			error: show404
		})
	}
}

// ------------------------------------------ //


function defaultsReplacement(cases, defaults) {
	let replacement = {}
	for (let [key, value] of Object.entries(defaults)) {
		replacement[key] = value
	}
	for (let [key, value] of Object.entries(cases)) {
		if (value !== undefined) {
			replacement[key] = value
		}
	}
	return replacement
}

function applyAttributes(object, cases, defaults) {
	for (let [key, value] of Object.entries(defaultsReplacement(cases, defaults))) {
		object.setAttribute(key, value)
	}
}

function generateAndApplyOptions(object, placeholder, cases) {
	placeholder = [[placeholder, '']]
	cases = placeholder.concat(Object.entries(cases))

	for (let [key, value] of cases) {
		let option = document.createElement('option')

		option.setAttribute('value', value)
		option.text = key

		object.appendChild(option)
	}
}

// ------------------------------------------ //


function generateTextInput(name, units, slots, data) {
	let wrap = document.createElement('div')
	wrap.className = 'form-group'

	if (String(slots.dementional) === 'true') {
		let span = document.createElement('div')
		span.className = 'form-group'

		let select = document.createElement('select')
		select.className = 'form-control'
		select.setAttribute('data-config-key', name + '.units')
		generateAndApplyOptions(select, 'выберите единицы измерения', units[slots.units])

		let _row = document.createElement('div')
		_row.className = 'row'
		let row = document.createElement('div')
		row.className = 'col-sm-12'

		for (let i = 0; i < 3; i++) {
			let col = document.createElement('div')
			col.className = 'col-sm-4'

			let input = document.createElement('input')
			applyAttributes(input, {'data-config-key': name + '.value[]'}, defaultTextInputCases)

			col.appendChild(input)
			row.appendChild(col)
		}

		span.appendChild(select)
		wrap.appendChild(span)
		_row.appendChild(row)
		wrap.appendChild(_row)
	} else {
		let textInput = document.createElement('input')
		applyAttributes(textInput, data, defaultTextInputCases)
		wrap.appendChild(textInput)
	}
	return wrap
}

function generateDropdown(name, slots) {
	let wrap = document.createElement('div')
	wrap.className = 'form-group'

	let dropdown = document.createElement('select')
	dropdown.className = 'form-control'
	dropdown.setAttribute('data-config-key', name)
	generateAndApplyOptions(dropdown, slots.title, slots.cases)

	wrap.appendChild(dropdown)
	return wrap
}

function generateCheckbox(name, slots) {
	let wrap = document.createElement('div')
	wrap.className = 'checkbox'
	let label = document.createElement('label')

	let checkbox = document.createElement('input')
	checkbox.type = 'checkbox'
	checkbox.setAttribute('data-config-key', name)
	checkbox.checked = (slots.default == 'True')

	label.appendChild(checkbox)
	label.appendChild(document.createTextNode(slots.title))
	wrap.appendChild(label)
	return wrap
}

function generateNumereticInput(name, units, slots) {
	let wrap = document.createElement('div')
	wrap.className = 'form-group'

	let label = document.createElement('label')
	label.appendChild(document.createTextNode(slots.title))
	wrap.appendChild(label)

	if (slots.dementional) {
		let span = document.createElement('div')
		span.className = 'form-group'

		let select = document.createElement('select')
		select.className = 'form-control'
		select.setAttribute('data-config-key', name + '.units')
		generateAndApplyOptions(select, 'выберите единицы измерения', units[slots.units])

		let _row = document.createElement('div')
		_row.className = 'row'
		let row = document.createElement('div')
		row.className = 'col-sm-12'

		for (let i = 0; i < 3; i++) {
			let col = document.createElement('div')
			col.className = 'col-sm-4'

			let input = document.createElement('input')
			applyAttributes(input, {
				'data-config-key': name + '.value[]',
				'minimum': slots.minimum,
				'maximum': slots.maximum,
				'step': slots.step,
			}, defaultNumereticInputCases)

			col.appendChild(input)
			row.appendChild(col)
		}

		span.appendChild(select)
		wrap.appendChild(span)
		_row.appendChild(row)
		wrap.appendChild(_row)
	} else {
		label.className = 'form-check-label'
		label.setAttribute('for', 'check-basic-collision')

		let mainInput = document.createElement('input')
		applyAttributes(mainInput, {
			'data-config-key': name + ((slots.units !== undefined) ? '.value' : ''),
			'minimum': slots.minimum,
			'maximum': slots.maximum,
			'step': slots.step,
		}, defaultNumereticInputCases)

		if (slots.units !== undefined) {
			let group = document.createElement('div')
			group.className = 'input-group'

			let span = document.createElement('span')
			let select = document.createElement('select')
			span.className = 'input-group-btn'
			select.className = 'btn btn-default'
			select.setAttribute('data-config-key', name + '.units')
			generateAndApplyOptions(select, 'выберите единицы измерения', units[slots.units])

			group.appendChild(mainInput)
			wrap.appendChild(group)

			span.appendChild(select)
			group.appendChild(span)
		} else {
			wrap.appendChild(mainInput)
		}
	}
	return wrap
}

// ------------------------------------------ //


function generateClassParameters(key, config) {
	let wrap = document.createElement('div')
	wrap.setAttribute('id', key)
	wrap.className = 'options-collapsable'

	for (let name of config.OBJECTS[key].cases) {
		let slots = config.CASES[name]
		if (slots.class === 'dropdown') {
			let dropdown = generateDropdown(name, slots)
			wrap.appendChild(dropdown)
		} else if (slots.class === 'checkbox') {
			let checkbox = generateCheckbox(name, slots)
			wrap.appendChild(checkbox)
		} else if (slots.class === 'number') {
			let numericInput = generateNumereticInput(name, config.UNITS, slots)
			wrap.appendChild(numericInput)
		} else if (slots.class === 'text') {
			let textInput = generateTextInput(name, config.UNITS, slots, {'data-config-key': 'PROBLEM'})
			wrap.appendChild(textInput)
		}
	}

	return wrap
}


function generateAllProblemForm(problemName, config) {
	let generalDiv = document.getElementsByClassName('row config-general')[0]

	let configNameField = generateTextInput(null, null, {}, {
		'data-config-key': 'PROBLEM',
		style: 'display: none',
		value: problemName
	})
	generalDiv.appendChild(configNameField)
	let general = config.GENERAL

	for (let name of general) {
		slots = config.CASES[name]
		name = 'GENERAL.' + name
		if (slots.class === 'dropdown') {
			let dropdown = generateDropdown(name, slots)
			generalDiv.appendChild(dropdown)
		} else if (slots.class === 'checkbox') {
			let checkbox = generateCheckbox(name, slots)
			generalDiv.appendChild(checkbox)
		} else if (slots.class === 'number') {
			let numericInput = generateNumereticInput(name, config.UNITS, slots)
			generalDiv.appendChild(numericInput)
		} else if (slots.class === 'text') {
			let textInput = generateTextInput(name, config.UNITS, slots, {'data-config-key': 'PROBLEM'})
			generalDiv.appendChild(textInput)
		}
	}

	if (config.OBJECTS !== undefined) {
		document.getElementById('config-objects').style = 'display: block;'
		let objectsSelect = document.getElementById('select-object-type')
		let objectsForm = document.getElementById('objects-form')

		for (let [key, value] of Object.entries(config.OBJECTS)) {
			let option = document.createElement('option')
			option.setAttribute('value', key)
			option.text = value.title
			objectsSelect.appendChild(option)
			objectsForm.append(generateClassParameters(key, config))
		}
	}

	autofocus()
}
