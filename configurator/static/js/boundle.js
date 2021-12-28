$.ajax({
    url: server_url + "/construct/configs/" + problem_name + ".yml",
    type: "GET",
    contentType: "application/yml; charset=utf-8",
    success: function (config) {
        generate_all(problem_name, jsyaml.load(config));
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
        console.log("poop!");
    });
});

let CONFIG = {};

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

function generate_select(cases) {
    let wrap = document.createElement("div");
    wrap.className = "cfg-hint-bx";

    let select = document.createElement("select");
    select.required = true;

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
    let basic = generate_basic_input(name, slots);

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

    let select = generate_select(slots.cases);
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

function generate_all(problem_name, config) {
    if (config.OBJECTS) {
        let save = document.getElementById("cfg-save");
        let next = document.getElementById("cfg-next");
        next.outerHTML = save.outerHTML;
        document.getElementById("cfg-page-3").remove();
    }
    let general = document.getElementById("cfg-general");

    for (let name of config.GENERAL) {
        let slots = config.CASES[name];
        name = "GENERAL." + name;

        if (slots.class === "dropdown") {
            let dropdown = generate_dropdown(name, slots);
            general.appendChild(dropdown);
        } else if (slots.class === "checkbox") {
            let checkbox = generate_checkbox(name, slots);
            general.appendChild(checkbox);
        } else if (slots.dementional) {
            let units = config.UNITS[slots.units];
            let input = generate_dimensional_input(name, slots, units);
            general.appendChild(input);
        } else if (slots.units) {
            let units = config.UNITS[slots.units];
            let input = generate_unit_input(name, slots, units);
            general.appendChild(input);
        } else {
            let input = generate_basic_input(name, slots);
            general.appendChild(input);
        }
    }
}
