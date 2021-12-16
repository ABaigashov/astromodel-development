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
