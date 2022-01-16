<style>
.cfg-check-item input:checked+span:after {
    background: url(<?php echo $server_url; ?>/construct/static/images/ch.svg) 50% no-repeat;
}
.cfg-select-inf-item .cfg-select:after {
    background: url(<?php echo $server_url; ?>/construct/static/images/arrow-down.svg) 50% no-repeat;
}
</style>

<link rel="preconnect" href="https://fonts.googleapis.com">
<link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
<link href="https://fonts.googleapis.com/css2?family=Open+Sans:wght@400;700&amp;display=swap" rel="stylesheet">
<script type="text/javascript" src="https://cdnjs.cloudflare.com/ajax/libs/jquery/3.5.0/jquery.min.js"></script>
<script type="text/javascript" src="https://cdnjs.cloudflare.com/ajax/libs/js-yaml/4.1.0/js-yaml.min.js"></script>
<script type="text/javascript" src="https://cdnjs.cloudflare.com/ajax/libs/jscolor/2.4.5/jscolor.min.js"></script>
<link rel="stylesheet" href="<?php echo $server_url; ?>/construct/static/css/config.css">
<script type="text/javascript" src="<?php echo $server_url; ?>/construct/static/js/configurator.js"></script>

<div class="cfg-container" id="cfg-page-1">
    <div class="cfg-table-bx">
        <div class="cfg-h3">Загрузка конфигурационного файла</div>
        <div class="cfg-step-bx">
            <div class="cfg-h4">Шаг 1. Создание пространства моделирования</div>
            <div class="cfg-file-bx">
                <div class="cfg-tt">
                    <div class="cfg-t1">Загрузите свой файл</div>
                    <div class="cfg-t2">Если у Вас уже есть конфигурациооный файл, то загрузите его (Форматы JSON, не
                        более 2 мб)</div>
                </div>
                <input type="file" id="cfg-file-input" style="display:none;" accept=".json" />
                <button class="cfg-btn" onclick="$('#cfg-file-input').trigger('click');">Прикрепить файл</button>
            </div>
            <div class="cfg-file-bx">
                <div class="cfg-tt">
                    <div class="cfg-t1">Создайте новый файл</div>
                    <div class="cfg-t2">Создайте новый конфигурациооный файл (Формат JSON)</div>
                </div>
                <a href="#cfg-page-2">
                    <button  class="cfg-btn">Создать файл</button>
                </a>
            </div>
            <div class="cfg-back-bx cfg-back-bx-cc">
                <button class="cfg-btn cfg-btn-framed" onclick="history.back();">Назад</button>
                <div class="cfg-tt">Если у Вас возникли сложности воспользуйтесь <a href="">инструкцией</a></div>
            </div>
        </div>
    </div>
</div>

<div class="cfg-container" id="cfg-page-2">
    <div class="cfg-table-bx">
        <div class="cfg-h3">Создание нового конфигурационного файла</div>
        <div class="cfg-h4">Шаг 1. Создание пространства моделирования</div>
        <div class="cfg-step-bx" id="cfg-general">

            <!-- general content here -->

        </div>
        <div class="cfg-back-bx">
            <a href="#cfg-page-1" class="cfg-btn cfg-btn-framed">Назад</a>
            <a href="#cfg-page-3" class="cfg-btn" id="cfg-next">Далее</a>
        </div>
    </div>
</div>


<div class="cfg-container" id="cfg-page-3">

    <div class="cfg-table-bx">
        <div class="cfg-h3">Создание объектов задачи</div>

        <div class="cfg-step-bx">

            <div class="cfg-table-wrap">
                <table class="cfg-table cfg-table-red">
                    <tbody id="cfg-table">
                        <tr>
                            <td>Название</td>
                            <td>Тип</td>
                            <td>Цвет</td>
                            <td></td>
                        </tr>
                    </tbody>
                </table>
            </div>
            <img class="cfg-plus-bt" src="<?php echo $server_url; ?>/construct/static/images/plus2.svg">

            <div class="cfg-back-bx">
                <a href="#cfg-page-2" class="cfg-btn cfg-btn-framed">Назад</a>
                <button class="cfg-btn" id="cfg-save">Cмоделировать</button>
            </div>
        </div>

    </div>

</div>


<div class="cfg-overlay"></div>

<div class="cfg-popup cfg-popup-form">
    <div class="cfg-close"><img src="<?php echo $server_url; ?>/construct/static/images/close.svg" alt=""></div>

    <div class="cfg-pop-form">
        <div class="cfg-h2">Настройки параметров объекта</div>
        <div id="cfg-objects-every">
            <div class="cfg-wr cfg-wr1">
                <div class="cfg-select-inf-item cfg-bx1">
                    <div class="cfg-lbl">
                        Название объекта
                        <div class="cfg-hint-bx">
                            <div class="cfg-opener"><img
                                    src="<?php echo $server_url; ?>/construct/static/images/info-circle.svg" alt=""></div>
                            <div class="cfg-frame">Пока ничего...</div>
                        </div>
                    </div>
                    <input type="text" placeholder="Текст" id="cfg-obj-name" data-config-key="name">
                </div>
                <div class="cfg-select-inf-item cfg-bx2 cfg-select-inf-item-color">
                    <div class="cfg-lbl">
                        Цвет объекта
                        <div class="cfg-hint-bx">
                            <div class="cfg-opener"><img
                                    src="<?php echo $server_url; ?>/construct/static/images/info-circle.svg" alt=""></div>
                            <div class="cfg-frame">Пока ничего...</div>
                        </div>
                    </div>
                    <input type="text" id="cfg-obj-color" data-config-key="color">
                </div>
            </div>
            <div class="cfg-select-inf-item">
                <div class="cfg-lbl">
                    Тип объекта
                    <div class="cfg-hint-bx">
                        <div class="cfg-opener"><img
                                src="<?php echo $server_url; ?>/construct/static/images/info-circle.svg" alt=""></div>
                        <div class="cfg-frame">Всплывающая информация по данному полю</div>
                    </div>
                </div>
                <div class="cfg-select">
                    <select id="cfg-object-type" required>
                        <option disabled selected hidden value="">Выберите</option>

                        <!-- objects types here -->

                    </select>
                </div>
            </div>
        </div>
        <div id="cfg-objects-content">

            <!-- objects content here -->

        </div>
        <div class="cfg-bt-wr">
            <button class="cfg-btn cfg-btn-framed">Удалить</button>
            <button class="cfg-btn" id="cfg-save-obj">Сохранить</Сохранить>
        </div>
    </div>
</div>
