var config = {
	objects: [
	],
};

const rgbToHex = (r, g, b) => '#' + [r, g, b]
  .map(x => x.toString(16).padStart(2, '0')).join('')

var list_button_template = '<button class="btn btn-xs btn-primary action-list" title="Список объектов" style="float: right;"><i class="fa fa-list"></i></button>';
var remove_button_template = '<button class="btn btn-xs btn-danger action-remove" title="Удалить объект" style="float: right;" data-id="{ID}"><i class="fa fa-minus"></i></button>';
var row_template = '<tr><td><button class="btn btn-xs btn-primary action-edit" title="Изменить настройки" data-id={ID}><i class="fa fa-edit"></i></button></td><td><div style="float: left; border: 1px solid #333; background: {COLOR};">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</div>&nbsp;&nbsp;{NAME}</td><td>{TYPE}</td><td>{REMOVE}</td></tr>';
var add_row_template = '<tr><td><button class="btn btn-xs btn-success action-add" title="Добавить объект"><i class="fa fa-plus"></i></button></td><td></td><td></td><td></td></tr>';
var editing_object_id = null;
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

		config.name = $('#config-name').val();
		config.dimensions = $('#config-dimensions').val();
		var blob = new Blob([JSON.stringify(config)], {type: "text/json;charset=utf-8"});
		//var FileSaver = new FileSaver();

		var config_file_name = $('#config-file').val();
		config_file_name = config_file_name.length > 0 ? config_file_name : 'config.json';
		config_file_name = config_file_name.endsWith('.json') ? config_file_name : config_file_name + '.json';
		window.saveAs(blob, config_file_name);
	});

	$('#file-import').on('change', function(e) {
		e.preventDefault();
		e.stopPropagation();
		ImportConfiguration(e, function(json){
			if(json) {
				config = json;
				$('#config-name').val(config.name ? config.name : '');
				$('#config-dimensions').val(config.dimensions ? config.dimensions : '');
				UpdateTable();
			}
		});
	});

	$('#btn-update').on('click', function(e) {
		e.preventDefault();
		e.stopPropagation();
		if(editing_object_id) {
			var new_object = ConfigureObject(editing_object_id);
			UpdateConfiguration(new_object);
			UpdateTable();
		}
		$('#modal-options').modal('hide');
	});

	UpdateTable();
});

/*Обновляет отображение попапа конфигурации объекта*/
function UpdateModalOptions(type) {
	var control_target = type + '-object';
	$('.options-collapsable').css('display', 'none');
	$('#' + control_target).css('display', 'block');
}

/*Обновляет отображение таблицы конфига*/
function UpdateTable() {
	var table_html = '';
	if(config.objects) {
		for(var s in config.objects) {
			var id = config.objects[s].id;
			var name = (config.objects[s].name) ? config.objects[s].name : '';
			var type = (config.objects[s].type) ? config.objects[s].type.description : '';
			var color = (config.objects[s].color) ? config.objects[s].color : [255,255,255];

			var remove_button = remove_button_template.replace('{ID}', id);

			var row_html = row_template;
			row_html = row_html.replace('{ID}', id);
			row_html = row_html.replace('{NAME}', name);
			row_html = row_html.replace('{TYPE}', type);
			row_html = row_html.replace('{REMOVE}', remove_button)
			row_html = row_html.replace('{COLOR}', rgbToHex(color[0], color[1], color[2]));

			table_html += row_html;
		}
	}

	table_html += add_row_template;
	$('#objects-list').html(table_html);

	$('.action-add').on('click', function(e){
		var object = {
			id: GetNewObjectId(3),
			name: 'Новый объект'
		}
		if(!config.objects) {
			config.objects = [];
		}
		config.objects.push(object);
		UpdateTable();
	});

	$('.action-list').on('click', function(e){

		$('#modal-list').modal();
	});

	$('.action-remove').on('click', function(e){
		var object_id = $(this).data('id');
		if(config.objects) {
			for(var s in config.objects) {
				if(config.objects[s].id === object_id) {
					console.log(s);
					config.objects.splice(s, 1);
					break;
				}
			}
		}
		UpdateTable();
	});

	$('.action-edit').on('click', function(e){
		var object_id = $(this).data('id');
		var object = GetObject(object_id);
		if(object != null) {
			editing_object_id = object_id;
			ApplyModalOptions(object);
			UpdateModalOptions(object.type ? object.type.value : '');
			$('#modal-options').modal();
		}
	});
}

/*Назначает значения характеристик объекта в поля ввода*/
function ApplyModalOptions(object) {
	$('#modal-options [data-config-key="name"]').val(object.name);
	$('#modal-options [data-config-key="type"]').val(object.type ? object.type.value : '');

	var color = object.color ? object.color : GetRandomColor();
	if(!jscolorPicker) {
		jscolorPicker = new jscolor($('#modal-options [data-config-key="color"]')[0]);
	}
	jscolorPicker.fromRGB(color[0], color[1], color[2]);

	if(object.quantity) {
		$('#modal-options [data-config-key="quantity"]').val(object.quantity[0]);
	}
	else {
		$('#modal-options [data-config-key="quantity"]').val('');
	}

	if(object.coords) {
		$('#modal-options [data-config-key="coord_x"]').val((object.coords.length > 0) ? object.coords[0] : '');
		$('#modal-options [data-config-key="coord_y"]').val((object.coords.length > 1) ? object.coords[1] : '');
		$('#modal-options [data-config-key="coord_z"]').val((object.coords.length > 2) ? object.coords[2] : '');
	}
	else {
		$('#modal-options [data-config-key="coord_x"]').val('');
		$('#modal-options [data-config-key="coord_y"]').val('');
		$('#modal-options [data-config-key="coord_z"]').val('');
	}

	if(object.speed) {
		$('#modal-options [data-config-key="speed_set"]').prop('checked', true);
		$('#modal-options [data-config-key="speed_x"]').val((object.speed.length > 0) ? object.speed[0] : '');
		$('#modal-options [data-config-key="speed_y"]').val((object.speed.length > 1) ? object.speed[1] : '');
		$('#modal-options [data-config-key="speed_z"]').val((object.speed.length > 2) ? object.speed[2] : '');
	}
	else {
		$('#modal-options [data-config-key="speed_set"]').prop('checked', false);
		$('#modal-options [data-config-key="speed_x"]').val('');
		$('#modal-options [data-config-key="speed_y"]').val('');
		$('#modal-options [data-config-key="speed_z"]').val('');
	}

	if(object.range) {
		$('#modal-options [data-config-key="range_min"]').val(object.range[0]);
		$('#modal-options [data-config-key="range_max"]').val(object.range[1]);
	}
	else {
		$('#modal-options [data-config-key="range_min"]').val('');
		$('#modal-options [data-config-key="range_max"]').val('');
	}

	if(object.electro) {
		$('#modal-options [data-config-key="charge_set"]').prop('checked', true);
		$('#modal-options [data-config-key="charge"]').val(object.electro.charge[0]);
		$('#modal-options [data-config-key="charge_min"]').val(object.electro.charge[0]);
		$('#modal-options [data-config-key="charge_max"]').val((object.electro.charge.length > 1) ? object.electro.charge[1] : object.electro.charge[0]);
	}
	else {
		$('#modal-options [data-config-key="charge_set"]').prop('checked', false);
		$('#modal-options [data-config-key="charge"]').val('');
		$('#modal-options [data-config-key="charge_min"]').val('');
		$('#modal-options [data-config-key="charge_max"]').val('');
	}

	if(object.gravity) {
		$('#modal-options [data-config-key="mass_set"]').prop('checked', true);
		$('#modal-options [data-config-key="mass"]').val(object.gravity.mass[0]);
		$('#modal-options [data-config-key="mass_min"]').val(object.gravity.mass[0]);
		$('#modal-options [data-config-key="mass_max"]').val((object.gravity.mass.length > 1) ? object.gravity.mass[1] : object.gravity.mass[0]);
	}
	else {
		$('#modal-options [data-config-key="mass_set"]').prop('checked', false);
		$('#modal-options [data-config-key="mass"]').val('');
		$('#modal-options [data-config-key="mass_min"]').val('');
		$('#modal-options [data-config-key="mass_max"]').val('');
	}

	if(object.collision) {
		$('#modal-options [data-config-key="radius_set"]').prop('checked', true);
		$('#modal-options [data-config-key="radius"]').val(object.collision.radius[0]);
		$('#modal-options [data-config-key="radius_min"]').val(object.collision.radius[0]);
		$('#modal-options [data-config-key="radius_max"]').val((object.collision.radius.length > 1) ? object.collision.radius[1] : object.collision.radius[0]);
	}
	else {
		$('#modal-options [data-config-key="radius_set"]').prop('checked', false);
		$('#modal-options [data-config-key="radius"]').val('');
		$('#modal-options [data-config-key="radius_min"]').val('');
		$('#modal-options [data-config-key="radius_max"]').val('');
	}
}

/*Создает новый объект из пользовательских данных*/
function ConfigureObject(object_id) {
	var object_type = $('#select-object-type').val();
	var target_section = object_type + '-object';
	var color = HexToRGB('#' + $('#modal-options [data-config-key="color"]').val());

	var object = {
		id: object_id,
		name: $('#modal-options [data-config-key="name"]').val(),
		color: color,
		type: {
			value: $('#select-object-type').val(),
			description: $('#select-object-type').find('option[value="' + object_type + '"]').data('description'),
		}
	}

	var coord_x = parseFloat($('#' + target_section + ' [data-config-key="coord_x"]').val());
	var coord_y = parseFloat($('#' + target_section + ' [data-config-key="coord_y"]').val());
	var coord_z = parseFloat($('#' + target_section + ' [data-config-key="coord_z"]').val());
	object.coords = [isNaN(coord_x) ? 0 : coord_x, isNaN(coord_y) ? 0 : coord_y, isNaN(coord_z) ? 0 : coord_z];

	var speed_x = parseFloat($('#' + target_section + ' [data-config-key="speed_x"]').val());
	var speed_y = parseFloat($('#' + target_section + ' [data-config-key="speed_y"]').val());
	var speed_z = parseFloat($('#' + target_section + ' [data-config-key="speed_z"]').val());
	object.speed = [isNaN(speed_x) ? 0 : speed_x, isNaN(speed_y) ? 0 : speed_y, isNaN(speed_z) ? 0 : speed_z];

	var charge = parseFloat($('#' + target_section + ' [data-config-key="charge"]').val());
	var val_min = parseFloat($('#' + target_section + ' [data-config-key="charge_min"]').val());
	var val_max = parseFloat($('#' + target_section + ' [data-config-key="charge_max"]').val());

	if(object_type === 'basic') {
		object.electro = {
			charge: [isNaN(charge) ? 0 : charge]
		};
	}
	else {
		val_min = isNaN(val_min) ? 0 : val_min;
		val_max = isNaN(val_max) ? 0 : val_max;
		object.electro = { charge: [Math.min(val_min, val_max), Math.max(val_min, val_max)] };
	}

	var mass = parseFloat($('#' + target_section + ' [data-config-key="mass"]').val());
	var val_min = parseFloat($('#' + target_section + ' [data-config-key="mass_min"]').val());
	var val_max = parseFloat($('#' + target_section + ' [data-config-key="mass_max"]').val());

	if(object_type === 'basic') {
		object.gravity = {
			mass: [isNaN(mass) ? 0 : mass]
		};
	}
	else {
		val_min = isNaN(val_min) ? 0 : val_min;
		val_max = isNaN(val_max) ? 0 : val_max;
		object.gravity = { mass: [Math.min(val_min, val_max), Math.max(val_min, val_max)] };
	}

	var radius = parseFloat($('#' + target_section + ' [data-config-key="radius"]').val());
	var val_min = parseFloat($('#' + target_section + ' [data-config-key="radius_min"]').val());
	var val_max = parseFloat($('#' + target_section + ' [data-config-key="radius_max"]').val());

	if(object_type === 'basic') {
		object.collision = {
			radius: [isNaN(radius) ? 0 : radius]
		};
	}
	else {
		val_min = isNaN(val_min) ? 0 : val_min;
		val_max = isNaN(val_max) ? 0 : val_max;
		object.collision = { radius: [Math.min(val_min, val_max), Math.max(val_min, val_max)] };
	}

	if(object_type !== 'basic') {
		var quantity = parseInt($('#' + target_section + ' [data-config-key="quantity"]').val());
		object.quantity = isNaN(quantity) ? 0 : quantity;

		if(object_type === 'array') {
			var val_min = parseFloat($('#' + target_section + ' [data-config-key="range_min"]').val());
			var val_max = parseFloat($('#' + target_section + ' [data-config-key="range_max"]').val());
			val_min = isNaN(val_min) ? 0 : val_min;
			val_max = isNaN(val_max) ? 0 : val_max;
			object.range = [Math.min(val_min, val_max), Math.max(val_min, val_max)];
		}
	}
	return object;
}

/* Сравнивает 2 объекта, возвращает True если объекты совпадают */
function CompareObjects(object_a, object_b) {
	if(!object_b.type)
		return false;

	if(object_a.type.value != object_b.type.value)
		return false;

	if(object_a.type.value !== 'basic') {
		if(object_a.quantity != object_b.quantity)
			return false;

		if(object_a.type.value === 'array') {
			if(object_a.range[0] != object_b.range[0] || object_a.range[1] != object_b.range[1])
				return false;
		}
	}
	return true;
}

/* Обновляет конфигурацию после внесения изменений пользователем */
function UpdateConfiguration(new_object) {
	var old_object = GetObject(new_object.id);
	if(old_object) {
		if(CompareObjects(new_object, old_object)) {
		}
		else {
		}

		for(var s in config.objects) {
			if(config.objects[s].id === new_object.id) {
				config.objects[s] = new_object;
				break;
			}
		}
	}
}

/* Возвращает объект по id*/
function GetObject(id) {
	if(config.objects) {
		for(var s in config.objects) {
			if(config.objects[s].id === id) {
				return config.objects[s];
			}
		}
	}
	return null;
}

function GetRangeValue(range) {
	if(Array.isArray(range)) {
		if(range.length == 0) {
			return 0;
		}
		else if(range.length == 1) {
			return range[0];
		}
		else {
			return range[0] + (range[1] - range[0]) * Math.random();
		}
	}
	else {
		return range;
	}
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
		if(config.objects) {
			for(var s in config.objects) {
				if(config.objects[s].id === id) {
					is_uniq = false;
					break;
				}
			}
		}
	} while(!is_uniq);
	return id;
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
