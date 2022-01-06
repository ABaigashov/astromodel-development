var STRUCTURE;

$.ajax({
    url: server_url + "/construct/configs/" + problem_name + ".yml",
    type: "GET",
    contentType: "application/yml; charset=utf-8",
    success: function (config) {
        STRUCTURE = jsyaml.load(config);
        generate_all();
    },
});

function create_new_object() {
    let table = document.getElementById("cfg-table");
    let row = document.createElement("tr");
    let count = table.children.length;
    row.className = "cfg-filled";
    for (let name of [
        "Солнце",
        "Горячее",
        "Звезда",
        "<div class='cfg-color cfg-color" + count + "'></div>",
        "<div class='cfg-red cfg-pop-form-opener'><img src='" +
            server_url +
            "/construct/static/images/red-bt.svg'></div>",
    ]) {
        let column = document.createElement("td");
        column.innerHTML = name;
        row.appendChild(column);
    }
    table.appendChild(row);
}

$(document).ready(function ($) {
    $(".cfg-pop-video").on("click", function (event) {
        event.preventDefault();
        $(".cfg-overlay, .cfg-popup-video").fadeIn();
    });
    $(".cfg-pop-form-opener").on("click", function (event) {
        event.preventDefault();
        $(".cfg-overlay, .cfg-popup-form").fadeIn();
    });
    $(".cfg-overlay, .cfg-popup .cfg-close").on("click", function (event) {
        $(".cfg-overlay, .cfg-popup").fadeOut();
    });
    $(".cfg-plus-bt").on("click", function (event) {
        event.preventDefault();
        create_new_object();
        $(".cfg-overlay, .cfg-popup-form").fadeIn();
        $(".cfg-pop-form-opener").on("click", function (event) {
            event.preventDefault();
            $(".cfg-overlay, .cfg-popup-form").fadeIn();
        });
    });
    $("#cfg-save").on("click", function (event) {
        saveConfig();
    });
});

let OBJECTS = {};

window.onload = () => {
    window.location.href = window.location.hash || "#cfg-page-1";
};

function generate_hint(slots) {
    let wrap = document.createElement("div");
    wrap.className = "cfg-hint-bx";

    let opener = document.createElement("div");
    opener.className = "cfg-opener";

    let opener_icon = document.createElement("img");
    opener_icon.src = server_url + "/construct/static/images/info-circle.svg";

    opener.appendChild(opener_icon);

    let frame = document.createElement("div");
    frame.className = "cfg-frame";
    frame.textContent = slots.hint || "Пока ничего...";

    wrap.appendChild(opener);
    wrap.appendChild(frame);

    return wrap;
}

function generate_label(slots) {
    let wrap = document.createElement("div");
    wrap.className = "cfg-lbl";

    let label = document.createTextNode(slots.title);
    wrap.appendChild(label);

    let hint = generate_hint(slots);
    wrap.append(hint);

    return wrap;
}

function generate_select(cases, name) {
    let wrap = document.createElement("div");
    wrap.className = "cfg-hint-bx";

    let select = document.createElement("select");
    select.required = true;
    select.setAttribute("data-config-key", name);

    let placeholder = document.createElement("option");
    placeholder.selected = true;
    placeholder.hidden = true;
    placeholder.disabled = true;
    placeholder.innerText = "Выберите";

    select.appendChild(placeholder);

    for (let [key, value] of Object.entries(cases)) {
        let option = document.createElement("option");

        option.value = value;
        option.innerText = key;

        select.appendChild(option);
    }
    wrap.appendChild(select);

    return wrap;
}

function generate_basic_input(name, slots) {
    let wrap = document.createElement("div");
    wrap.className = "cfg-select-inf-item";

    let label = generate_label(slots);
    wrap.appendChild(label);

    let input = document.createElement("input");
    input.type = "text";
    input.placeholder = "Введите";
    input.setAttribute("data-config-key", name);

    wrap.appendChild(input);

    return wrap;
}

function generate_unit_input(name, slots, units) {
    let basic = generate_basic_input(name + ".value", slots);

    let dropdown = generate_dropdown(name + ".units", {
        title: "Ед. изм.",
        cases: units,
        hint: undefined,
    });
    dropdown.className = "cfg-sel-over";

    basic.appendChild(dropdown);

    return basic;
}

function generate_dimensional_input(name, slots, units) {
    let basic = generate_basic_input(name, slots);

    let dropdown = generate_dropdown(name + ".units", slots, units);
    dropdown.className = "cfg-sel-over";

    basic.appendChild(dropdown);

    return basic;
}

function generate_dropdown(name, slots) {
    let wrap = document.createElement("div");
    wrap.className = "cfg-select-inf-item";

    let label = generate_label(slots);
    wrap.appendChild(label);

    let select = generate_select(slots.cases, name);
    select.className = "cfg-select";
    wrap.appendChild(select);

    return wrap;
}

function generate_checkbox(name, slots) {
    let wrap = document.createElement("div");
    wrap.className = "cfg-ch-group";

    let label = document.createElement("label");
    label.className = "cfg-check-item";

    let checkbox = document.createElement("input");
    checkbox.type = "checkbox";
    checkbox.setAttribute("data-config-key", name);
    checkbox.checked = slots.default == "True";

    let span = document.createElement("span");
    span.textContent = slots.title;

    label.appendChild(checkbox);
    label.appendChild(span);

    let hint = generate_hint(slots);

    wrap.appendChild(label);
    wrap.appendChild(hint);

    return wrap;
}

function generate_all() {
    let general = document.getElementById("cfg-general");

    for (let name of STRUCTURE.GENERAL) {
        let slots = STRUCTURE.CASES[name];
        name = "GENERAL." + name;

        if (slots.class === "dropdown") {
            let dropdown = generate_dropdown(name, slots);
            general.appendChild(dropdown);
        } else if (slots.class === "checkbox") {
            let checkbox = generate_checkbox(name, slots);
            general.appendChild(checkbox);
        } else if (slots.dementional) {
            let units = STRUCTURE.UNITS[slots.units];
            let input = generate_dimensional_input(name, slots, units);
            general.appendChild(input);
        } else if (slots.units) {
            let units = STRUCTURE.UNITS[slots.units];
            let input = generate_unit_input(name, slots, units);
            general.appendChild(input);
        } else {
            let input = generate_basic_input(name, slots);
            general.appendChild(input);
        }
    }

    if (!STRUCTURE.OBJECTS) {
        let save = document.getElementById("cfg-save");
        let next = document.getElementById("cfg-next");
        next.outerHTML = save.outerHTML;
        document.getElementById("cfg-page-3").remove();
        return;
    }

    let objects_select = document.getElementById("cfg-object-type");
    console.log(objects_select);
}

function myencode(data) {
    data = JSON.stringify(data);
    data = encodeURI(data);
    data = btoa(data);
    return data;
}

function saveConfig() {
    // let cfg = {};
    // for (let name of STRUCTURE.GENERAL) {
    //     cfg.OBJECTS = OBJECTS;
    // }

    let cfg = JSON.parse(
        '{"OBJECTS":{"point_objects":[{"id":"xqy-fyl","index":2,"coords":{"units":"149597870700","value":["1","0","0"]},"speed":{"units":"1000","value":["0","30","0"]},"mass":{"units":"1","value":"10"},"charge":{"units":"","value":"0"},"radius":{"units":"1000","value":"700000"},"delay":{"units":"60","value":"0"},"trajectory":"","name":"Earth","color":[159,215,223],"type":{"value":"point_objects","description":"Единичный объект"}},{"id":"wxk-xfo","index":1,"coords":{"units":"1000","value":["0","0","0"]},"speed":{"units":"0.001","value":["0","0","0"]},"mass":{"units":"1.989e+30","value":"1"},"charge":{"units":"1","value":"0"},"radius":{"units":"1000","value":"700000"},"delay":{"units":"1","value":"0"},"trajectory":"","name":"Sun","color":[178,26,128],"type":{"value":"point_objects","description":"Единичный объект"}}],"EM_fields":[],"G_fields":[],"random_generators":[{"id":"uiu-pjr","index":3,"radius_scale":{"units":"1000","value":"7000000"},"mass_scale":{"units":"1000","value":"10"},"charge_scale":{"units":"0.001","value":"10"},"delay_scale":{"units":"1","value":"0"},"coordinate_scale":{"units":"149597870700","value":"1"},"velocity_scale":{"units":"1000","value":"20"},"particals_number":"5","equal_parameters":"","name":"gen1","color":[87,11,171],"type":{"value":"random_generators","description":"Генератор произвольного распределения объектов"}},{"id":"acj-jrs","index":4,"radius_scale":{"units":"1000","value":"9000000"},"mass_scale":{"units":"0.001","value":"1"},"charge_scale":{"units":"0.001","value":"23"},"delay_scale":{"units":"0.001","value":"0"},"coordinate_scale":{"units":"149597870700","value":"1"},"velocity_scale":{"units":"1000","value":"32"},"particals_number":"5","name":"gen2","color":[23,150,11],"type":{"value":"random_generators","description":"Генератор произвольного распределения объектов"}}],"ring_generator":[{"id":"vbz-lmc","index":5,"radius_scale":{"units":"1000","value":"9000000"},"mass_scale":{"units":"0.001","value":"1"},"charge_scale":{"units":"1","value":"1"},"delay_scale":{"units":"86400","value":"0"},"mass_center":{"units":"1.989e+30","value":"1"},"coords_center":{"units":"1000","value":["0","0","0"]},"speed_center":{"units":"0.000001","value":["0","0","0"]},"ring_radius":{"units":"149597870700","value":"1"},"point_velocity":{"units":"1000","value":"30"},"particals_number":"20","equal_parameters":"","name":"ring1","color":[194,0,0],"type":{"value":"ring_generator","description":"Генератор кругового распределения объектов"}}],"ellipse_generators":[{"id":"hkn-hoe","index":6,"radius_scale":{"units":"1000","value":"9000000"},"mass_scale":{"units":"0.001","value":"1"},"charge_scale":{"units":"0.001","value":"2"},"delay_scale":{"units":"0.000001","value":"0"},"mass_center":{"units":"1.989e+30","value":"1"},"coords_center":{"units":"149597870700","value":["0","0","0"]},"speed_center":{"units":"1000","value":["0","0","0"]},"major_axis":{"units":"149597870700","value":"2"},"focus":"-1","apside_angle":{"units":"0.017453292519943295","value":"70"},"particals_number":"30","exentricity":"0.8","equal_parameters":"","name":"ell1","color":[213,40,155],"type":{"value":"ellipse_generators","description":"Генератор эллиптического распределения объектов"}}]},"GENERAL":{"dimensions":"2","K":"1","output_graphics":"matplotlib","step":{"units":"86400","value":"1"},"edge":{"units":"149597870700","value":"2"},"name":"test","steps_number":"130","fps":"30","frames_gap":"1","scale_faktor":"1","gravity_point_interaction":"on","trajectory":"on"},"PROBLEM":"particle_simulator"}'
    );

    console.log(myencode(cfg));
    alert("saved");
}
