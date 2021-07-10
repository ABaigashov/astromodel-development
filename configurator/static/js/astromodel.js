var config = {
	OBJECTS: {
		point_objects: [],
	}
};

const rgbToHex = (r, g, b) => '#' + [r, g, b]
  .map(x => x.toString(16).padStart(2, '0')).join('')

var list_button_template = '<button class="btn btn-xs btn-primary action-list" title="Список объектов" style="float: right;"><i class="fa fa-list"></i></button>';
var remove_button_template = '<button class="btn btn-xs btn-danger action-remove" title="Удалить объект" style="float: right;" data-id="{ID}"><i class="fa fa-minus"></i></button>';
var row_template = '<tr><td><button class="btn btn-xs btn-primary action-edit" title="Изменить настройки" data-id={ID}><i class="fa fa-edit"></i></button></td><td><div style="float: left; border: 1px solid #333; background: {COLOR};">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</div>&nbsp;&nbsp;{NAME}</td><td>{TYPE}</td><td>{REMOVE}</td></tr>';
var add_row_template = '<tr><td><button class="btn btn-xs btn-success action-add" title="Добавить объект"><i class="fa fa-plus"></i></button></td><td></td><td></td><td></td></tr>';
var editing_object_id = null;
var editing_object_index = 0;
var jscolorPicker = null;

$(function(){
	$('.options-collapsable').css('display', 'none');

	$('#select-object-type').on('change', function(e){
		var value_selected = $(this).val();
		UpdateModalOptions(value_selected);
	});

	$('#btn-export').on('click', function(e) {
		e.preventDefault();
		e.stopPropagation();

		var new_config = {
			OBJECTS: {}
		};
		ApplyInputData(new_config, $('.config-general'));
		var objects_lists = GetListsOfObjects();
		for(var i in objects_lists) {
			var list_name = objects_lists[i];
			if(config.OBJECTS[list_name]) {
				new_config.OBJECTS[list_name] = config.OBJECTS[list_name];
			}
		}
		config = new_config;

		var blob = new Blob([JSON.stringify(config)], {type: "text/json;charset=utf-8"});
		//var FileSaver = new FileSaver();

		var config_file_name = $('#config-file').val();
		config_file_name = config_file_name.length > 0 ? config_file_name : 'config.json';
		config_file_name = config_file_name.endsWith('.json') ? config_file_name : config_file_name + '.json';
		config_file_name = config_file_name.replace(/\ /g, '_');
		window.saveAs(blob, Translit(config_file_name));
	});

	$('#file-import').on('change', function(e) {
		e.preventDefault();
		e.stopPropagation();
		ImportConfiguration(e, function(json){
			if(json) {
				config = json;
				RestoreConfigObjects();
				FillInputForm(config, $('.config-general'));
				UpdateTable();
			}
		});
	});

	$('#btn-update').on('click', function(e) {
		e.preventDefault();
		e.stopPropagation();
		if(editing_object_id) {
			var new_object = ConfigureObject(editing_object_id, editing_object_index);
			UpdateConfiguration(new_object);
			UpdateTable();
		}
		$('#modal-options').modal('hide');
	});

	UpdateTable();
});


/*Обновляет отображение попапа конфигурации объекта*/
function UpdateModalOptions(type) {
	$('.options-collapsable').css('display', 'none');
	if(type != '') {
		$('#' + type).css('display', 'block');
	}
}

/*Обновляет отображение таблицы конфига*/
function UpdateTable() {
	var table_html = '';
	var list = GetObjectsList();

	for(var s in list) {
		var id = list[s].id;
		var name = (list[s].name) ? list[s].name : '';
		var type = (list[s].type) ? list[s].type.description : '';
		var color = (list[s].color) ? list[s].color : [255,255,255];

		var remove_button = remove_button_template.replace('{ID}', id);

		var row_html = row_template;
		row_html = row_html.replace('{ID}', id);
		row_html = row_html.replace('{NAME}', name);
		row_html = row_html.replace('{TYPE}', type);
		row_html = row_html.replace('{REMOVE}', remove_button)
		row_html = row_html.replace('{COLOR}', rgbToHex(color[0], color[1], color[2]));

		table_html += row_html;
	}

	table_html += add_row_template;
	$('#objects-list').html(table_html);

	$('.action-add').on('click', function(e){
		var object = {
			id: GetNewObjectId(3),
			name: 'Новый объект',
			index: GetNewObjectIndex()
		}

		var objects_lists = GetListsOfObjects();
		var config_objects = config.OBJECTS;
		for(var i in objects_lists) {
			var list_name = objects_lists[i];
			if(!config_objects[list_name]) {
				config_objects[list_name] = [];
			}
		}

		config_objects.point_objects.push(object);
		config.OBJECTS = config_objects;
		UpdateTable();
	});

	$('.action-list').on('click', function(e){
		$('#modal-list').modal();
	});

	$('.action-remove').on('click', function(e){
		var object_id = $(this).data('id');
		var objects_lists = GetListsOfObjects();
		var config_objects = config.OBJECTS;
		for(var i in objects_lists) {
			var list_name = objects_lists[i];
			if(config_objects[list_name]) {
				for(var s in config.OBJECTS[list_name]) {
					if(config_objects[list_name][s].id === object_id) {
						config_objects[list_name].splice(s, 1);
						break;
					}
				}
			}
		}
		config.OBJECTS = config_objects;
		UpdateTable();
	});

	$('.action-edit').on('click', function(e){
		var object_id = $(this).data('id');
		var object = GetObject(object_id);
		if(object != null) {
			editing_object_id = object_id;
			editing_object_index = object.index;
			if(!editing_object_index) {
				editing_object_index = GetNewObjectIndex();
			}
			ApplyModalOptions(object);
			UpdateModalOptions(object.type ? object.type.value : '');
			$('#modal-options').modal();
		}
	});
}
/*Writes into struct input data from source element */
function ApplyInputData(struct, inputs_root) {
	ApplyInputDataList(struct, $(inputs_root).find('select'));
	ApplyInputDataList(struct, $(inputs_root).find('input[type="text"]'));
	ApplyInputDataList(struct, $(inputs_root).find('input[type="number"]'));
	ApplyInputDataList(struct, $(inputs_root).find('input[type="checkbox"]').filter((index, item) => {
			return $(item).prop('checked');
		})
	);
}
/*Writes into struct input data from list of input elements*/
function ApplyInputDataList(struct, inputs) {
	inputs.each((index, item) => {
		UpdateStructure(struct, item);
	});
}
/*Writes into struct input data from specific input element*/
function UpdateStructure(struct, elem) {
	var name = $(elem).data('config-key');
	var value = $(elem).val();
	UpdateStructureFiled(struct, name, value);
}

/*Writes into structvalue by specific path*/
function UpdateStructureFiled(struct, path, value) {
	var level = path.split('.');
	var name = level[0];

	if(name.endsWith('[]')) {
		name = name.substr(0, name.length - 2);
		if(!(name in struct)) {
			struct[name] = [];
		}
		struct[name].push(value);
	}
	else {
		if(level.length > 1) {
			if(!(name in struct)) {
				struct[name] = {};
			}
			level.splice(0, 1);
			UpdateStructureFiled(struct[name], level.join('.'), value);
		}
		else {
			struct[name] = value;
		}
	}
}

/*Clear form */
function ClearInputForm(inputs_root) {
	$(inputs_root).find('input').val('');
	$(inputs_root).find('select').val('');
	$(inputs_root).find('input').prop('checked', false);
}
/*Fills struct into form */
function FillInputForm(struct, inputs_root) {
	var fieldList = GetConfigValuesList(struct);
	for(var key in fieldList) {
		var value = fieldList[key];
		if(Array.isArray(value)) {
			$(inputs_root).find('[data-config-key="' + key + '"]').each((index, item) => {
				if($(item).attr('type') == 'checkbox') {
					var inp_val = $(item).attr('value');
					for(var i in value) {
						if(value[i] == inp_val) {
							$(item).prop('checked', true);
							break;
						}
					}
				}
				else {
					$(item).val(value[index]);
				}
			});

			// Fix of assing of array to single fields inside each collapsalbe options set.
			// Otherwise first value will be assigned to first input and ingnored in anoher options sets.
			$(inputs_root).find('.options-collapsable').each((index, options_set) => {
				$(options_set).find('[data-config-key="' + key + '"]').each((index, item) => {
					if($(item).attr('type') == 'checkbox') {
						var inp_val = $(item).attr('value');
						for(var i in value) {
							if(value[i] == inp_val) {
								$(item).prop('checked', true);
								break;
							}
						}
					}
					else {
						$(item).val(value[index]);
					}
				});
			});
		}
		else {
			if($(inputs_root).find('[data-config-key="' + key + '"]').attr('type') == 'checkbox') {
				$(inputs_root).find('[data-config-key="' + key + '"]').prop('checked', true);
			}
			else {
				$(inputs_root).find('[data-config-key="' + key + '"]').val(value);
			}
		}
	}
}

/*Restores list of objects for config file of old version*/
function RestoreConfigObjects() {
	if(!config.OBJECTS) {
		var config_objects = {};
		$(GetListsOfObjects()).each(function(index, item){
			if(config[item]) {
				config_objects[item] = config[item];
			}
		});
		config.OBJECTS = config_objects;
	}
}

/*Return list of fileld of provided level of the config file*/
function GetConfigValuesList(struct, prefix) {
	prefix = (prefix) ? prefix : '';
	var list = {};
	for(var key in struct) {
		var value = struct[key];
		if(Array.isArray(value)) {
			list[prefix + key + '[]'] = value;
		}
		else if(value instanceof Object) {
			list = {...list, ...GetConfigValuesList(value, prefix + key + '.')};
		}
		else {
			list[prefix + key] = value;
		}
	}
	return list;
}

/*Назначает значения характеристик объекта в поля ввода*/
function ApplyModalOptions(object) {
	ClearInputForm($('#modal-options'));
	FillInputForm(object, $('#modal-options'));
	var color = object.color ? object.color : GetRandomColor();
	if(!jscolorPicker) {
		jscolorPicker = new jscolor($('#modal-options [data-config-key="color"]')[0]);
	}
	jscolorPicker.fromRGB(color[0], color[1], color[2]);
}

/*Создает новый объект из пользовательских данных*/
function ConfigureObject(object_id, object_index) {
	var new_object = {
		id: object_id,
		index: object_index
	}
	var object_type = $('#select-object-type').val();
	ApplyInputData(new_object, $('#' + object_type));

	var color = HexToRGB('#' + $('#modal-options [data-config-key="color"]').val());
	new_object.name = $('#modal-options [data-config-key="name"]').val();
	new_object.color = color;
	new_object.type = {
		value: object_type,
		description: $('#select-object-type').find('option[value="' + object_type + '"]').html(),
	};
	return new_object;
}

/* Обновляет конфигурацию после внесения изменений пользователем */
function UpdateConfiguration(new_object) {
	var objects_lists = GetListsOfObjects();
	var config_objects = config.OBJECTS;
	for(var i in objects_lists) {
		var list_name = objects_lists[i];
		if(config_objects[list_name]) {
			for(var s in config_objects[list_name]) {
				if(config_objects[list_name][s].id === new_object.id) {
					config_objects[list_name].splice(s, 1);
					break;
				}
			}
		}
	}

	if(!config_objects[new_object.type.value]) {
		config_objects[new_object.type.value] = [];
	}
	config_objects[new_object.type.value].push(new_object);
	config.OBJECTS = config_objects;
}

/* Возвращает объект по id*/
function GetObject(id) {
	var objects_lists = GetListsOfObjects();
	var config_objects = config.OBJECTS;
	for(var i in objects_lists) {
		var list_name = objects_lists[i];
		if(config_objects[list_name]) {
			for(var s in config_objects[list_name]) {
				if(config_objects[list_name][s].id === id) {
					return config_objects[list_name][s];
				}
			}
		}
	}
	return null;
}

function GetListsOfObjects() {
	var list = [];
	$('#select-object-type option').each((index, item) => {
		var list_name = $(item).attr('value');
		if(list_name) {
			if(list.indexOf(list_name) < 0)
				list.push(list_name);
		}
	});
	return list;
}

function GetObjectsList() {
	var list = [];
	var config_objects = config.OBJECTS;
	$(GetListsOfObjects()).each(function(index, item){
		if(config_objects[item]) {
			list = list.concat(config_objects[item]);
		}
	});
	list.sort(function(a,b){
		return a.index - b.index;
	});
	return list;
}

function ImportConfiguration(e, callback) {
    var files = e.target.files;
    for (var i = 0, f; f = files[i]; i++) {
		if (!f.type.match('application/json')) {
			continue;
		}
		$('#config-file').val(files[i].name);

		var reader = new FileReader();
		reader.onload = (function(theFile) {
			return function(evt) {
				var result = evt.target.result;
				if(result.startsWith('data:application/json;base64')) {
					result = result.substr(29);
					result =  window.Base64.decode(result);
				}

				var json_config = null;
				try {
					json_config = JSON.parse(result);
				}
				catch(err) {
					console.log(err);
				}
				callback(json_config);
			};
		})(f);

		// Read in the image file as a data URL.
		reader.readAsDataURL(f);
	}
}

/*Возвращает уникальный идентификатор в рамках текущего конфига*/
function GetNewObjectId(size) {
	var id = '';
	var is_uniq = false;
	do {
		id = GenRandomId(size);
		is_uniq = true;
		var objects_lists = GetListsOfObjects();
		for(var i in objects_lists) {
			var list_name = objects_lists[i];
			if(config[list_name]) {
				for(var s in config[list_name]) {
					if(config[list_name][s].id === id) {
						is_uniq = false;
						break;
					}
				}
			}
		}
	} while(!is_uniq);
	return id;
}
/*Возвращает уникальный идентификатор в рамках текущего конфига*/
function GetNewObjectIndex(size) {
	var object_index = 0;
	var list = GetObjectsList();
	for(var i in list) {
		var item = list[i];
		if(item.index) {
			object_index = Math.max(object_index, item.index);
		}
	}
	return object_index + 1;
}
/*Возвращает случайный идентификатор*/
function GenRandomId(size) {
	var characterSet = [
		'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j',
		'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 't',
		'u', 'v', 'w', 'x', 'y', 'z'
	];
	var id = '';
	var index = 0;
	for(var i = 0; i < size; i++) {
		index = Math.floor(Math.random() * characterSet.length);
		id += characterSet[index];
	}
	id += '-';
	for(var i = 0; i < size; i++) {
		index = Math.floor(Math.random() * characterSet.length);
		id += characterSet[index];
	}
	return id;
}

function GetRandomColor() {
	var color = [Math.floor(Math.random() * 255), Math.floor(Math.random() * 255), Math.floor(Math.random() * 255)];
	console.log('Random color: ' + color);
	return color;
}

function HexToRGB(hex){
    var c;
    if(/^#([A-Fa-f0-9]{3}){1,2}$/.test(hex)){
        c= hex.substring(1).split('');
        if(c.length== 3){
            c= [c[0], c[0], c[1], c[1], c[2], c[2]];
        }
        c= '0x'+c.join('');
        return [(c>>16)&255, (c>>8)&255, c&255];
    }
    throw new Error('Bad Hex');
}

function RGBToHex(color) {
	var c = [Math.floor(color[0] * 255), Math.floor(color[1] * 255), Math.floor(color[2] * 255)];
	return ((1 << 24) + (c[0] << 16) + (c[1] << 8) + c[2]).toString(16).padStart(2, '0').join('');
}

function Translit(str) {
    var dict = {
        'а': 'a', 'б': 'b', 'в': 'v', 'г': 'g', 'д': 'd',
        'е': 'e', 'ё': 'e', 'ж': 'j', 'з': 'z', 'и': 'i',
        'к': 'k', 'л': 'l', 'м': 'm', 'н': 'n', 'о': 'o',
        'п': 'p', 'р': 'r', 'с': 's', 'т': 't', 'у': 'u',
        'ф': 'f', 'х': 'h', 'ц': 'c', 'ч': 'ch', 'ш': 'sh',
        'щ': 'shch', 'ы': 'y', 'э': 'e', 'ю': 'u', 'я': 'ya'
    }, n_str = [];

    str = str.replace(/[ъь]+/g, '').replace(/й/g, 'i');

    for ( var i = 0; i < str.length; ++i ) {
       n_str.push(
              dict[ str[i] ]
           || dict[ str[i].toLowerCase() ] == undefined && str[i]
           || dict[ str[i].toLowerCase() ].toUpperCase()
       );
    }

    return n_str.join('');
}
