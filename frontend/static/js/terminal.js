$(document).ready(($) => {
    $("#cfg-file-input").on("change", (event) => {
        event.preventDefault();
        let file = $("#cfg-file-input").prop("files")[0];
        let reader = new FileReader();
        reader.readAsText(file, "UTF-8");
        reader.onload = (event) => {
            api_request(
                "job_execute",
                {
                    config: JSON.parse(event.target.result),
                    filename: file.name.split(".")[0],
                    session,
                },
                () => {}
            );
        };
    });
});

const session = Array(64)
    .fill()
    .map((x) => "0123456789abcdef"[(16 * Math.random()) | 0])
    .join("");

function api_request(url, data, callback) {
    let date = new Date();
    $.ajax({
        method: "POST",
        url: request_server + "/api/" + url,
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
                callback(data.answer);
            }
        })
        .fail(function (jqXHR, textStatus) {
            console.log(jqXHR);
            alert(textStatus);
        });
}

function update_list() {
    api_request("job_list", { session: session }, (answer) => {
        let rows =
            "<tr><td>№</td><td>Сценарий</td><td>Тип проблемы</td><td>Состояние</td><td>Создана</td><td>Завершена</td><td></td></tr>";
        for (let i in answer) {
            let row = "<tr class='cfg-filled'>";
            row += "<td>" + answer[i].uid + "</td>";
            row += "<td>" + answer[i].name + "</td>";
            row += "<td>" + answer[i].problem + "</td>";

            if (answer[i].progress == -1) {
                let user = $(this).data("user");
                let job = $(this).data("job");
                let resultUrl =
                    request_server +
                    "/static/job/" +
                    answer[i].user_uid +
                    "/" +
                    answer[i].uid +
                    "/error.txt";
                row +=
                    '<td><a class="cfg-pop-pop" href="' +
                    encodeURI(resultUrl) +
                    '" style="color:#ff0000;">Показать ошибку</a></td>';
            } else if (answer[i].progress == 1) {
                let user = $(this).data("user");
                let job = $(this).data("job");
                let resultUrl =
                    request_server +
                    "/static/job/" +
                    answer[i].user_uid +
                    "/" +
                    answer[i].uid +
                    "/result" +
                    answer[i].extention;
                row +=
                    '<td><a class="cfg-pop-pop" href="' +
                    encodeURI(resultUrl) +
                    '" style="color:#0000ff;">Показать результат</a></td>';
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
                '<td><a href="" class="cfg-del" data-uid="' +
                answer[i].uid +
                '"><img src="' +
                server_url +
                '/construct/static/images/del.svg" alt=""></a></td>';
            row += "</tr>";
            rows += row;
        }

        $("#cfg-table").html(rows);
        if ($(".cfg-container").css("opacity") === "0") {
            $(".cfg-container").css("opacity", "1");
        }

        setTimeout(() => {
            update_list();
        }, 1000);

        $(".cfg-pop-pop").on("click", async (event) => {
            event.preventDefault();

            $(".cfg-overlay, .cfg-popup-video").fadeIn();

            let responce = event.target.href;

            let rows =
                '<div class="cfg-close"><img src="' +
                server_url +
                '/construct/static/images/close.svg" alt=""></div><div class="cfg-vid-bx">';

            let ext = responce.split(".").slice(-1)[0];
            let error = false;

            switch (ext) {
                case "mp4":
                    rows +=
                        '<video autoplay loop style="width:100%;max-height:65vh"><source src="' +
                        responce +
                        '" type="video/mp4"></video>';
                    break;
                case "txt":
                    try {
                        let text = await $.ajax({
                            method: "GET",
                            url: responce,
                        });
                        rows += "<pre>" + text + "</pre>";
                    } catch {
                        error = true;
                        rows +=
                            '<pre style="line-height:20;text-align:center;margin:50px">Файл с логом ошибки отсутствует.</pre>';
                    }
                    break;
                default:
                    rows +=
                        '<pre style="line-height:20;text-align:center;margin:50px">Неизвестный формат данных ".' +
                        ext +
                        '"</pre>';
            }

            rows += '</div><div class="cfg-btn-wr">';
            rows +=
                '<a href="' +
                responce
                    .split("/")
                    .slice(0, -1)
                    .concat("config.json")
                    .join("/") +
                '" class="cfg-btn cfg-btn-framed cfg-save-config">Сохранить конфигурацию</a>';

            if (error) {
                rows +=
                    '<a href="' +
                    responce
                        .split("/")
                        .slice(0, -1)
                        .concat("config.json")
                        .join("/") +
                    '" class="cfg-btn cfg-btn-filled cfg-retry-config">Редактировать конфигурацию</a>';
            } else {
                rows +=
                    '<a href="' +
                    responce
                        .split("/")
                        .slice(0, -1)
                        .concat("config.json")
                        .join("/") +
                    '" class="cfg-btn cfg-btn-framed cfg-retry-config">Редактировать конфигурацию</a>';
                if (event.target.innerText === "Показать результат") {
                    rows +=
                        '<a href="' +
                        responce +
                        '" class="cfg-btn cfg-btn-filled cfg-save-result">Сохранить результат</a>';
                } else {
                    let = rows +=
                        '<a href="' +
                        responce +
                        '" class="cfg-btn cfg-btn-filled cfg-save-result">Сохранить ошибку</a>';
                }
            }
            rows += "</div>";

            $("#cfg-result").html(rows);

            $(".cfg-popup .cfg-close").on("click", (event) => {
                event.preventDefault();
                $(".cfg-overlay, .cfg-popup").fadeOut();
            });

            $(".cfg-retry-config").on("click", (event) => {
                event.preventDefault();
                fetch(event.target.href)
                    .then((res) => res.json())
                    .then((json) => {
                        window.name = btoa(encodeURI(JSON.stringify(json)));
                        setTimeout(() => {
                            window.location.pathname = "/" + json.PROBLEM + "/";
                        }, 100);
                    });
            });

            $(".cfg-save-result, .cfg-save-config").on(
                "click",
                async (event) => {
                    event.preventDefault();
                    let resource = event.target.href;

                    let type = await fetch(resource, {
                        method: "HEAD",
                    }).then((res) => {
                        return res.headers.get("content-type");
                    });

                    let fix_utf = ["application/json", "plain/text"];

                    if (fix_utf.includes(type)) {
                        fetch(resource)
                            .then((res) => res.text())
                            .then((text) => {
                                let blob = new Blob([eval("`" + text + "`")], {
                                    type,
                                });
                                window.saveAs(
                                    blob,
                                    resource.split("/").slice(-1)[0]
                                );
                            });
                    } else {
                        fetch(resource)
                            .then((res) => res.blob())
                            .then((blob) => {
                                window.saveAs(
                                    blob,
                                    resource.split("/").slice(-1)[0]
                                );
                            });
                    }
                }
            );
        });

        $(".cfg-del img").on("click", (event) => {
            event.preventDefault();
            event.stopPropagation();

            let uid =
                event.target.parentElement.parentElement.parentElement
                    .children[0].innerText;

            if (
                confirm(
                    "Вы уверены что хотите навсегда удалить задачу №" +
                        uid +
                        " из вашего списка?"
                )
            ) {
                api_request(
                    "job_remove",
                    {
                        job_uniq: uid,
                        session: session,
                    },
                    () => null
                );
            }
        });
    });
}

function formatDate(timestamp, format) {
    if (timestamp == 0) {
        return "Никогда";
    }
    let datetime = new Date();
    datetime.setTime(timestamp * 1000);

    let date = [];
    date.push(datetime.getFullYear());
    date.push(datetime.getMonth() + 1);
    date.push(datetime.getDate());
    date.push(datetime.getHours());
    date.push(datetime.getMinutes());
    date.push(datetime.getSeconds());

    for (let i in date) {
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

try {
    try {
        request_server;
    } catch (e) {
        let article = document.getElementsByClassName("article")[0];
        article.outerHTML = `<div style="background:white;width:100%;height:60vh;display:flex;flex-direction:column;justify-content:center;align-items:center;font-family:'Ubuntu Mono',monospace;"><h1>Error -1</h1><h3>No server url</h3></div>`;
        throw Error("");
    }

    fetch(request_server)
        .then((result) => {
            if (!result.ok) {
                throw Error("");
            }

            let login = document
                .getElementsByClassName("user-snippet")[0]
                .getAttribute("data-user-email");

            api_request("user_register", { login, name: "" }, () => {
                api_request("user_authorize", { login, session }, async () => {
                    if (window.name) {
                        let name = window.name;
                        try {
                            let cfg = JSON.parse(decodeURI(atob(name)));
                            let date = new Date();
                            let data = (
                                await $.ajax({
                                    method: "POST",
                                    url: request_server + "/api/job_list",
                                    data: JSON.stringify({
                                        timezone: date.getTimezoneOffset(),
                                        data: { session },
                                    }),
                                    dataType: "json",
                                })
                            ).answer.map((x) => x.name);

                            setTimeout(() => {
                                let filename;

                                if (cfg.GENERAL.name) {
                                    filename = cfg.GENERAL.name;
                                } else {
                                    filename = 1;
                                    while (
                                        data.includes("Задача_" + filename)
                                    ) {
                                        filename += 1;
                                    }
                                    filename = "Задача_" + filename;
                                }

                                api_request(
                                    "job_execute",
                                    {
                                        config: cfg,
                                        filename,
                                        session,
                                    },
                                    () => {}
                                );
                                console.log("Config loaded via configurator");
                                window.name = "";
                            }, 100);
                        } catch (e) {
                            console.log(e);
                        }
                    }
                    update_list();
                });
            });
        })
        .catch(() => {
            let article = document.getElementsByClassName("article")[0];
            article.outerHTML = `<div style="background:white;width:100%;height:60vh;display:flex;flex-direction:column;justify-content:center;align-items:center;font-family:'Ubuntu Mono',monospace;"><h1>Error -1</h1><h3>Wrong server url</h3></div>`;
        });
} catch {}
