var STRUCTURE;
var jscolorPicker;

$.ajax({
    url: server_url + "/construct/configs/" + problem_name + ".yml",
    type: "GET",
    contentType: "application/yml; charset=utf-8",
    success: function (config) {
        STRUCTURE = jsyaml.load(config);
        generate_all();
    },
});

$(document).ready(($) => {
    $(".cfg-pop-form-opener").on("click", (event) => {
        event.preventDefault();
        $(".cfg-overlay, .cfg-popup-form").fadeIn();
    });
    $(".cfg-popup .cfg-close").on("click", (event) => {
        event.preventDefault();
        let children = [...document.getElementById("cfg-table").children].slice(
            1
        );
        for (let child of children) {
            let inner = JSON.parse(child.children[1].innerText);
            if (inner.edit && !inner.type) {
                child.remove();
            } else {
                delete inner.edit;
                child.children[1].innerText = JSON.stringify(inner);
            }
        }
        $(".cfg-overlay, .cfg-popup").fadeOut();
    });
    $(".cfg-popup .cfg-btn-framed").on("click", (event) => {
        event.preventDefault();
        let children = [...document.getElementById("cfg-table").children].slice(
            1
        );
        for (let child of children) {
            if (JSON.parse(child.children[1].innerText).edit) {
                child.remove();
            }
        }
        $(".cfg-overlay, .cfg-popup").fadeOut();
    });
    $(".cfg-plus-bt").on("click", (event) => {
        event.preventDefault();
        new_obj();
    });
    $("#cfg-save").on("click", (event) => {
        event.preventDefault();
        save_config();
    });
    $("#cfg-model").on("click", (event) => {
        event.preventDefault();
        model_config();
    });
    $("#cfg-save-obj").on("click", (event) => {
        event.preventDefault();
        save_obj();
    });
    $("#cfg-file-input").on("change", (event) => {
        event.preventDefault();
        let file = $("#cfg-file-input").prop("files")[0];
        let reader = new FileReader();
        reader.readAsText(file, "UTF-8");
        reader.onload = (event) => {
            load_config(JSON.parse(event.target.result));
            setTimeout(() => window.scrollTo(0, 0), 1);
        };
    });
});

let OBJECTS = {};

window.onload = () => {
    if (
        ["", "#", "#cfg-page-1", "#cfg-page-2", "#cfg-page-3"].includes(
            window.location.hash
        )
    ) {
        window.location.href = window.location.hash || "#cfg-page-1";
    }
    if (window.name) {
        let name = window.name;
        setTimeout(() => {
            try {
                let cfg = JSON.parse(decodeURI(atob(name)));
                try {
                    load_config(cfg);
                    console.log("Config loaded via terminal");
                } catch (e) {
                    console.log(cfg);
                    console.error(e);
                }
                window.name = "";
                setTimeout(() => window.scrollTo(0, 0), 1);
            } catch {}
        }, 100);
    }
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

    if (slots.short) {
        wrap.classList.add("cfg-lbl-short");
    }

    let label = document.createTextNode(slots.title);
    wrap.appendChild(label);

    let hint = generate_hint(slots);
    wrap.append(hint);

    return wrap;
}

function generate_select(cases, name, _default) {
    let wrap = document.createElement("div");
    wrap.className = "cfg-hint-bx";

    let select = document.createElement("select");
    select.required = true;
    select.setAttribute("data-config-key", name);

    let placeholder = document.createElement("option");

    if (_default === undefined) {
        placeholder.selected = true;
        placeholder.value = "";
        placeholder.innerText = "Выберите";
    } else {
        placeholder.selected = true;
        placeholder.value = cases[_default];
        placeholder.innerText = _default;
    }

    select.appendChild(placeholder);
    for (let [key, value] of Object.entries(cases)) {
        if (_default !== undefined && key !== _default) {
            let option = document.createElement("option");

            option.value = value;
            option.innerText = key;

            select.appendChild(option);
        }
    }
    console.log(select);
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
    slots.short = true;
    let basic = generate_basic_input(name + ".value", slots);

    let input = basic.getElementsByTagName("input")[0];

    input.type = "number";
    input.placeholder = "Число";

    let dropdown = generate_dropdown(
        name + ".units",
        {
            title: "Ед. изм.",
            cases: units.values,
            hint: units.hint,
        },
        units.default
    );
    dropdown.className = "cfg-sel-over";

    basic.appendChild(dropdown);

    return basic;
}

function generate_dimensional_input(name, slots, units) {
    let wrap = document.createElement("div");
    wrap.className = "cfg-select-inf-item";

    let label = generate_label(slots);
    wrap.appendChild(label);

    let dim = document.createElement("div");
    dim.className = "cfg-dim-bx";

    for (let i = 0; i < 3; i++) {
        let box = document.createElement("div");
        box.className = "cfg-bx";

        let letter = document.createElement("div");
        letter.className = "cfg-let";
        letter.innerText = String.fromCharCode(88 + i);

        box.appendChild(letter);

        let input = document.createElement("input");
        input.type = "number";
        input.placeholder = "Число";
        input.setAttribute("data-config-key", name + ".value[]");

        box.appendChild(input);

        dim.appendChild(box);
    }
    wrap.appendChild(dim);

    let dropdown = generate_dropdown(
        name + ".units",
        {
            title: "Ед. изм.",
            cases: units.values,
            hint: units.hint,
        },
        units.default
    );
    dropdown.className = "cfg-line-bx";

    wrap.appendChild(dropdown);

    return wrap;
}

function generate_dropdown(name, slots, _default) {
    let wrap = document.createElement("div");
    wrap.className = "cfg-select-inf-item";

    let label = generate_label(slots);
    wrap.appendChild(label);

    let select = generate_select(slots.cases, name, _default);
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

function general_field(wrap, name, slots) {
    if (slots.class === "dropdown") {
        let dropdown = generate_dropdown(name, slots);
        wrap.appendChild(dropdown);
    } else if (slots.class === "checkbox") {
        let checkbox = generate_checkbox(name, slots);
        wrap.appendChild(checkbox);
    } else if (slots.dementional) {
        let units = STRUCTURE.UNITS[slots.units];
        let input = generate_dimensional_input(name, slots, units);
        wrap.appendChild(input);
    } else if (slots.units) {
        let units = STRUCTURE.UNITS[slots.units];
        let input = generate_unit_input(name, slots, units);
        wrap.appendChild(input);
    } else {
        let input = generate_basic_input(name, slots);
        wrap.appendChild(input);
    }
}

function generate_object(wrap, cases) {
    for (let name of cases) {
        general_field(wrap, name, STRUCTURE.CASES[name]);
    }
}

function generate_all() {
    let general = document.getElementById("cfg-general");

    for (let name of STRUCTURE.GENERAL) {
        general_field(general, "GENERAL." + name, STRUCTURE.CASES[name]);
    }

    if (!STRUCTURE.OBJECTS) {
        let save = document.getElementById("cfg-save");
        let next = document.getElementById("cfg-next");
        next.outerHTML = save.outerHTML;
        document.getElementById("cfg-page-3").remove();
        return;
    }

    let types = document.getElementById("cfg-object-type");
    for (let name of Object.keys(STRUCTURE.OBJECTS)) {
        let type = document.createElement("option");

        type.value = name;
        type.innerText = STRUCTURE.OBJECTS[name].title;

        types.appendChild(type);
    }

    let objects = document.getElementById("cfg-objects-content");

    for (let [name, values] of Object.entries(STRUCTURE.OBJECTS)) {
        let content = document.createElement("div");

        content.id = "cfg-obj-" + name;
        generate_object(content, values.cases);
        content.style.display = "none";

        objects.appendChild(content);
    }

    $("#cfg-object-type").change((event) => {
        for (let div of objects.children) {
            div.style.display = "none";
        }
        document.getElementById("cfg-obj-" + event.target.value).style.display =
            "block";
    });

    let color_input = $("#cfg-obj-color")[0];
    jscolor.presets.default = {
        previewSize: Number(getComputedStyle(color_input).height.slice(0, -2)),
    };
    jscolorPicker = new jscolor(color_input);
}

function encodeSelector(el) {
    let path = [],
        parent;
    while ((parent = el.parentNode)) {
        path.unshift(
            `${el.tagName}:${[].indexOf.call(parent.children, el) + 1}`
        );
        el = parent;
    }
    return `${path.join(">")}`.toLowerCase();
}

function decodeSelector(code) {
    code = code.replaceAll(":", ":nth-child(").replaceAll(">", ") > ") + ")";
    return document.querySelector(code);
}

function get_all_dataed(element) {
    let data = [];
    for (let part of element.children) {
        if (part.hasAttribute("data-config-key")) {
            data.push(part);
        } else if (part.childElementCount !== 0) {
            data = data.concat(get_all_dataed(part));
        }
    }
    return data;
}

function gen_data_by_element(element, exclude) {
    let data = {};
    let selectors = [];

    let elements = get_all_dataed(element);
    if (exclude === undefined) {
        elements = elements.concat(
            get_all_dataed(document.getElementById("cfg-objects-every"))
        );
    }

    for (let element of elements) {
        let name = element.getAttribute("data-config-key").split(".");
        for (let i = 0; i < name.length - 1; i++) {
            let ex = "data." + name.slice(0, i + 1).join(".");
            eval(ex + " = " + ex + " || {}");
        }
        let value = element.value;
        if (element.type === "checkbox") {
            selectors.push(
                "decodeSelector('" +
                    encodeSelector(element) +
                    "').checked = " +
                    element.checked
            );
            value = element.checked;
        } else if (element.id === "cfg-obj-color") {
            selectors.push("jscolorPicker.fromString('" + value + "')");
            value = parseInt(value.slice(1), 16);
            value = [(value >> 16) & 255, (value >> 8) & 255, value & 255];
        } else {
            if (element.tagName === "SELECT" && element.value === "Выберите") {
                selectors.push(
                    "decodeSelector('" +
                        encodeSelector(element) +
                        "').selectedIndex = 0"
                );
            } else {
                selectors.push(
                    "decodeSelector('" +
                        encodeSelector(element) +
                        "').value = '" +
                        value +
                        "'"
                );
            }
            if (element.type === "number" || name.slice(-1)[0] === "units") {
                value = Number(value);
            }
        }

        if (name.slice(-1)[0].endsWith("[]")) {
            let ex = "data." + name.join(".").slice(0, -2);
            eval(ex + " = (" + ex + " || []).concat(value)");
        } else {
            eval("data." + name.join(".") + " = value");
        }
    }

    return [selectors, data];
}

function clear_inputs() {
    let form = document.getElementsByClassName("cfg-pop-form")[0];
    for (let input of form.getElementsByTagName("input")) {
        if (input.type === "checkbox") {
            input.checked = false;
        } else {
            input.value = "";
        }
    }
    for (let input of form.getElementsByTagName("select")) {
        input.selectedIndex = 0;
    }
}

function update_hooks() {
    let objects = document.getElementById("cfg-objects-content");

    $(".cfg-pop-form-opener").on("click", function (event) {
        event.preventDefault();
        let row = event.target.parentElement.parentElement.parentElement;
        let data = JSON.parse(row.children[1].innerText);
        let selectors = JSON.parse(row.children[0].innerText);

        data.edit = true;
        row.children[1].innerText = JSON.stringify(data);
        $(".cfg-overlay, .cfg-popup-form").fadeIn();

        clear_inputs();

        document.getElementById("cfg-object-type").value = data.type;

        for (let div of objects.children) {
            div.style.display = "none";
        }
        document.getElementById("cfg-obj-" + data.type).style.display = "block";

        for (let selector of selectors) {
            eval(selector);
        }
    });
}

function new_obj() {
    $(".cfg-overlay, .cfg-popup-form").fadeIn();

    clear_inputs();

    let table = document.getElementById("cfg-table");
    let name = document.getElementById("cfg-obj-name");

    let count = table.children.length;
    name.value = "Объект_" + count;

    let type = document.getElementById("cfg-object-type");
    type.selectedIndex = 0;

    let objects = document.getElementById("cfg-objects-content");
    for (let div of objects.children) {
        div.style.display = "none";
    }
    let color = "#" + (((1 << 24) * Math.random()) | 0).toString(16);
    jscolorPicker.fromString(color);

    let row = document.createElement("tr");
    row.className = "cfg-filled";

    let id = GenRandomId(3);

    let cases = Object.entries([
        "",
        `{"edit":true, "id":"${id}"}`,
        name.value,
        "",
        "<div class='cfg-color' style='background-color:" + color + "'></div>",
        "<div class='cfg-red cfg-pop-form-opener'><img src='" +
            server_url +
            "/construct/static/images/red-bt.svg'></div>",
    ]);

    for (let [index, name] of cases) {
        let column = document.createElement("td");
        if (index < 2) {
            column.style.display = "none";
        }
        column.innerHTML = name;
        row.appendChild(column);
    }
    table.appendChild(row);

    update_hooks();
}

function save_obj() {
    let types = document.getElementById("cfg-object-type");
    if (types.selectedIndex === 0) {
        alert("Выберите тип объекта");
        return;
    }
    let [selectors, data] = gen_data_by_element(
        document.getElementById("cfg-obj-" + types.value)
    );

    $(".cfg-overlay, .cfg-popup-form").fadeOut();

    let table = document.getElementById("cfg-table");
    let children = [...table.children].slice(1);

    for (let child of children) {
        let inner = JSON.parse(child.children[1].innerText);
        if (inner.edit) {
            delete inner.edit;

            data.type = types.value;
            child.children[0].innerText = JSON.stringify(selectors);
            child.children[1].innerText = JSON.stringify({
                ...data,
                id: inner.id,
            });

            child.children[2].innerText = data.name;
            child.children[3].innerText = data.type;
            child.children[4].innerHTML =
                "<div class='cfg-color' style='background-color:rgb(" +
                data.color.join(", ") +
                ")'></div>";
        }
    }
}

function GenRandomId(size) {
    var characterSet = [
        "a",
        "b",
        "c",
        "d",
        "e",
        "f",
        "g",
        "h",
        "i",
        "j",
        "k",
        "l",
        "m",
        "n",
        "o",
        "p",
        "q",
        "r",
        "s",
        "t",
        "u",
        "v",
        "w",
        "x",
        "y",
        "z",
    ];
    var id = "";
    var index = 0;
    for (var i = 0; i < size; i++) {
        index = Math.floor(Math.random() * characterSet.length);
        id += characterSet[index];
    }
    id += "-";
    for (var i = 0; i < size; i++) {
        index = Math.floor(Math.random() * characterSet.length);
        id += characterSet[index];
    }
    return id;
}

function load_config(data) {
    for (let selector of data.SELECTORS[0]) {
        eval(selector);
    }
    let selectors = {};

    for (let [type, objects] of Object.entries(data.OBJECTS)) {
        for (let obj of objects) {
            let index = obj.index;
            delete obj.index;
            obj.type = type;
            selectors[index] = {
                selectors: data.SELECTORS[index + 1],
                data: obj,
            };
        }
    }

    for (let obj of Object.values(selectors)) {
        let table = document.getElementById("cfg-table");
        let row = document.createElement("tr");
        row.className = "cfg-filled";

        let color =
            "#" +
            (
                (1 << 24) +
                (obj.data.color[0] << 16) +
                (obj.data.color[1] << 8) +
                obj.data.color[2]
            )
                .toString(16)
                .slice(1);

        let cases = Object.entries([
            JSON.stringify(obj.selectors),
            JSON.stringify(obj.data),
            obj.data.name,
            obj.data.type,
            "<div class='cfg-color' style='background-color:" +
                color +
                "'></div>",
            "<div class='cfg-red cfg-pop-form-opener'><img src='" +
                server_url +
                "/construct/static/images/red-bt.svg'></div>",
        ]);

        for (let [index, name] of cases) {
            let column = document.createElement("td");
            if (index < 2) {
                column.style.display = "none";
            }
            column.innerHTML = name;
            row.appendChild(column);
        }
        table.appendChild(row);
    }

    update_hooks();
    window.location.href = "#cfg-page-2";
}

function __get_config() {
    let general = document.getElementById("cfg-general");
    let [selectors, config] = gen_data_by_element(general, true);

    config.SELECTORS = [selectors];
    config.PROBLEM = problem_name;
    config.OBJECTS = {};
    for (let type of Object.keys(STRUCTURE.OBJECTS)) {
        config.OBJECTS[type] = [];
    }

    let table = document.getElementById("cfg-table");
    let children = [...table.children].slice(1);

    for (let [index, child] of Object.entries(children)) {
        data = JSON.parse(child.children[1].innerText);
        selectors = JSON.parse(child.children[0].innerText);
        let type = data.type;
        delete data.type;
        data.index = Number(index);
        config.SELECTORS.push(selectors);
        config.OBJECTS[type].push(data);
    }
    config.CFG_REF = window.location.href;
    return config;
}

function save_config() {
    let config = __get_config();

    let blob = new Blob([JSON.stringify(config)], {
        type: "application/json",
    });
    let filename = "config";
    if (config.GENERAL.name) {
        filename = cfg.GENERAL.name;
    }
    window.saveAs(blob, filename + ".json");
}

function model_config() {
    let config = __get_config();

    let enc = btoa(encodeURI(JSON.stringify(config)));
    window.name = enc;
    setTimeout(() => {
        window.location.href = window.location.origin + terminal_url;
    }, 10);
}
