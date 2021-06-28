
function autofocus() {
	let general = [...document.getElementsByClassName('row config-general')[0].children]
	general = [general[1]].concat(general.slice(3))
	let elements = []
	general.forEach(div => {
		if (div.children[0].className == 'form-check-label') {
			if (div.children[1].className == 'input-group') {
				elements.push(div.children[1].children[0], div.children[1].children[1].children[0])
			} else {
				elements.push(div.children[1])
			}
		} else {
			elements.push(div.children[0])
		}
	})
	elements[0].focus()
	for (let i = 0; i < elements.length - 1; i++) {
		elements[i].addEventListener('keydown', function(event) {
			if (event.keyCode === 13) {
				elements[i + 1].focus()
			}
		})
	}
	elements.slice(-1)[0].addEventListener('keydown', function(event) {
		if (event.keyCode === 13) {
			document.activeElement.blur()
		}
	})
}