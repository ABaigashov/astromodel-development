$.ajax({
    url: server_url + "/construct/configs/" + problem_name + ".yml",
    type: "GET",
    contentType: "application/yml; charset=utf-8",
    success: function (config) {
        generate_all(problem_name, jsyaml.load(config));
    },
});

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
});

window.onload = () => {
    window.location.href = window.location.hash || "#cfg-page-1";
};

function generate_all(problem_name, config) {
    console.log(problem_name);
    console.log(config);
}

function generate_hint(slots) {
    let wrap = document.createElement("div");
    wrap.className = "cfg-hint-bx";

    let opener = document.createElement("div");
    opener.className = "cfg-opener";

    let opener_icon = document.createElement("img");
    opener_icon.src =
        "<?php echo $server_url; ?>/construct/static/images/info-circle.svg";

    opener.appendChild(opener_icon);

    let frame = document.createElement("div");
    frame.className = "cfg-frame";
    frame.textContent = slots.hint || "Пока ничего...";

    wrap.appendChild(opener);
    wrap.appendChild(frame);

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

/*
<div class="cfg-select-inf-item">
    <div class="cfg-lbl">
        Название сценария
        <div class="cfg-hint-bx">
            <div class="cfg-opener"><img src="https://off.ddns.net:8008/construct/static/images/info-circle.svg" alt="">
            </div>
            <div class="cfg-frame">Всплывающая информация по данному полю</div>
        </div>
    </div>
    <input type="text" placeholder="Текст">
</div>
*/
