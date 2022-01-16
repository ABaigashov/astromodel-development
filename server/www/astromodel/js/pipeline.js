var step = 1;
var userSession = "";
var userConfig = null;
//var apiUrl = 'http://127.0.0.1:8080/';
// var apiUrl = 'https://off.ddns.net:8888/';
var apiUrl = window.location.origin + "/";

window.addEventListener("message", function (event) {
    if (
        event.origin != "https://astromodel.ru" &&
        !event.origin.startsWith("http://localhost")
    ) {
        alert("wrong domain: " + event.origin);
        return;
    }

    if (event.data.command == "auth") {
        authorize(event.data.user, event.data.name);
    }
});

$(function () {
    setStage(step);

    $("#step-1 form").on("submit", function (e) {
        e.preventDefault();
        e.stopPropagation();
        var login = $(this).find('input[name="login"]').val();
        var name = $(this).find('input[name="name"]').val();
        var data = {
            login: login,
            name: name,
        };
        apiRequest(apiUrl + "api/", "user_register", data, function (answer) {
            userSession = randomUniq(64);
            var data = {
                login: login,
                session: userSession,
            };
            apiRequest(
                apiUrl + "api/",
                "user_authorize",
                data,
                function (answer) {
                    setStage(step + 1);
                    updateJobList();
                }
            );
        });
    });
    $('#step-2 form input[type="file"]').on("change", function (e) {
        e.preventDefault();
        e.stopPropagation();
        importConfiguration(e, function (json) {
            if (json) {
                userConfig = json;
                var fileName = $('#step-2 form input[type="file"]').val();
                fileName = fileName.replace(/\\/g, "/");
                fileName = fileName.substr(fileName.lastIndexOf("/") + 1);
                fileName = fileName.substr(0, fileName.lastIndexOf("."));

                var data = {
                    session: userSession,
                    config: userConfig,
                    filename: fileName,
                };

                $("#astro-pending").css("display", "block");
                apiRequest(
                    apiUrl + "api/",
                    "job_execute",
                    data,
                    function (answer) {
                        $('#step-2 form input[type="file"]').val("");
                        $("#astro-pending").css("display", "none");
                    }
                );
            }
        });
    });
});

function authorize(login, name) {
    var data = {
        login: login,
        name: name,
    };

    $("#astro-pending").css("display", "block");
    apiRequest(apiUrl + "api/", "user_register", data, function (answer) {
        userSession = randomUniq(64);
        var data = {
            login: login,
            session: userSession,
        };
        apiRequest(apiUrl + "api/", "user_authorize", data, function (answer) {
            $("#astro-pending").css("display", "none");
            setStage(step + 1);
            updateJobList();
        });
    });
}

function updateJobList() {
    var data = {
        session: userSession,
    };
    apiRequest(apiUrl + "api/", "job_list", data, function (answer) {
        var rows = "";
        for (var i in answer) {
            var row = "<tr>";
            row += "<td>" + answer[i].uid + "</td>";
            row += "<td>" + answer[i].name + "</td>";
            row += "<td>" + answer[i].problem + "</td>";

            if (answer[i].progress == -1) {
                var user = $(this).data("user");
                var job = $(this).data("job");
                var resultUrl =
                    apiUrl +
                    "static/job/" +
                    answer[i].user_uid +
                    "/" +
                    answer[i].uid +
                    "/error.txt";
                row +=
                    '<td><a href="' +
                    resultUrl +
                    '" target="_blank" class="get-result" style="color: #da343f;">Показать ошибку</a></td>';
            } else if (answer[i].progress == 1) {
                var user = $(this).data("user");
                var job = $(this).data("job");
                var resultUrl =
                    apiUrl +
                    "static/job/" +
                    answer[i].user_uid +
                    "/" +
                    answer[i].uid +
                    "/result" +
                    answer[i].extention;
                row +=
                    '<td><a href="' +
                    resultUrl +
                    '" target="_blank" class="get-result">Показать результат</a></td>';
            } else if (answer[i].progress > 0) {
                row +=
                    "<td>" +
                    parseFloat((answer[i].progress * 100).toFixed(2)) +
                    "%</td>";
            } else {
                switch (answer[i].state) {
                    case "await":
                        row += "<td>В очереди</td>";
                        break;
                    default:
                        row += "<td>" + answer[i].state + "</td>";
                        break;
                }
            }

            row +=
                "<td>" +
                formatDate(answer[i].date_created, "Y/M/D H:I:S") +
                "</td>";
            row +=
                "<td>" +
                formatDate(answer[i].date_finished, "Y/M/D H:I:S") +
                "</td>";
            row +=
                '<td><button class="btn btn-sm btn-danger btn-remove" data-uid="' +
                answer[i].uid +
                '"><i class="fa fa-times"></i></button></td>';
            row += "</tr>";
            rows += row;
        }
        $("#job-list").html(rows);

        var tm = setTimeout(function () {
            updateJobList();
        }, 1000);

        $(".btn-remove").on("click", function (e) {
            e.preventDefault();
            e.stopPropagation();

            clearTimeout(tm);
            var uid = $(this).data("uid");
            if (
                confirm(
                    "Вы уверены что хотите навсегда удалить задачу #" +
                        uid +
                        " из вашего списка?"
                )
            ) {
                var data = {
                    job_uniq: uid,
                    session: userSession,
                };
                $("#astro-pending").css("display", "block");
                apiRequest(
                    apiUrl + "api/",
                    "job_remove",
                    data,
                    function (answer) {
                        $("#astro-pending").css("display", "none");
                        updateJobList();
                    }
                );
            } else {
                updateJobList();
            }
        });
    });
}

function setStage(stage) {
    $("#step-1, #step-2, #step-3, #step-4, #step-5").css("display", "none");
    $("#step-" + stage).css("display", "block");
    step = stage;
}

function formatDate(timestamp, format) {
    if (timestamp == 0) {
        return "Никогда";
    }
    var datetime = new Date();
    datetime.setTime(timestamp * 1000);

    var date = [];
    date.push(datetime.getFullYear());
    date.push(datetime.getMonth() + 1);
    date.push(datetime.getDate());
    date.push(datetime.getHours());
    date.push(datetime.getMinutes());
    date.push(datetime.getSeconds());

    for (var i in date) {
        date[i] = date[i] < 10 ? "0" + date[i] : date[i].toString();
    }

    format = format.replace(/Y/g, date[0]);
    format = format.replace(/M/g, date[1]);
    format = format.replace(/D/g, date[2]);
    format = format.replace(/H/g, date[3]);
    format = format.replace(/I/g, date[4]);
    format = format.replace(/S/g, date[5]);
    return format;
}

function apiRequest(url, method, data, done) {
    var date = new Date();
    $.ajax({
        method: "POST",
        url: url + method,
        data: JSON.stringify({
            data: data,
            timezone: date.getTimezoneOffset(),
        }),
        dataType: "json",
    })
        .done(function (data) {
            if (data.redirect) {
                window.location.href = data.redirect;
            } else if (data.error) {
                alert(data.error);
            } else {
                done(data.answer);
            }
        })
        .fail(function (jqXHR, textStatus) {
            console.log(jqXHR);
            alert(textStatus);
        });
}

function importConfiguration(e, callback) {
    var files = e.target.files;
    for (var i = 0, f; (f = files[i]); i++) {
        if (!f.type.match("application/json")) {
            continue;
        }

        var reader = new FileReader();
        reader.onload = (function (theFile) {
            return function (evt) {
                var result = evt.target.result;
                if (result.startsWith("data:application/json;base64")) {
                    result = result.substr(29);
                    result = window.Base64.decode(result);
                }

                var json_config = null;
                try {
                    json_config = JSON.parse(result);
                } catch (err) {
                    console.log(err);
                }
                callback(json_config);
            };
        })(f);

        // Read in the image file as a data URL.
        reader.readAsDataURL(f);
    }
}

function randomUniq(length) {
    var chars = [
        "0",
        "1",
        "2",
        "3",
        "4",
        "5",
        "6",
        "7",
        "8",
        "9",
        "a",
        "b",
        "c",
        "d",
        "e",
        "f",
    ];
    var uniq = "";
    while (uniq.length < length) {
        var c = Math.floor(Math.random() * chars.length);
        uniq += chars[c];
    }
    return uniq;
}
