<style>
.cfg-check-item input:checked+span:after {
    background: url(<?php echo $server_url; ?>/construct/static/images/ch.svg) 50% no-repeat;
}

.cfg-select-inf-item .cfg-select:after {
    background: url(<?php echo $server_url; ?>/construct/static/images/arrow-down.svg) 50% no-repeat;
}
</style>

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
                <button class="cfg-btn">Прикрепить файл</button>
            </div>
            <div class="cfg-file-bx">
                <div class="cfg-tt">
                    <div class="cfg-t1">Создайте новый файл</div>
                    <div class="cfg-t2">Создайте новый конфигурациооный файл (Формат JSON)</div>
                </div>
                <a href="#cfg-page-2" class="cfg-btn">Создать файл</a>
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

        <div class="cfg-step-bx">
            <div class="cfg-h4">Шаг 1. Создание пространства моделирования</div>

            <div class="cfg-select-inf-item">
                <div class="cfg-lbl">
                    Название сценария
                    <div class="cfg-hint-bx">
                        <div class="cfg-opener"><img
                                src="<?php echo $server_url; ?>/construct/static/images/info-circle.svg" alt="">
                        </div>
                        <div class="cfg-frame">Всплывающая информация по данному полю</div>
                    </div>
                </div>
                <input type="text" placeholder="Текст">
            </div>
            <div class="cfg-select-inf-item">
                <div class="cfg-lbl">
                    Размерность пространства
                    <div class="cfg-hint-bx">
                        <div class="cfg-opener"><img
                                src="<?php echo $server_url; ?>/construct/static/images/info-circle.svg" alt="">
                        </div>
                        <div class="cfg-frame">Всплывающая информация по данному полю</div>
                    </div>
                </div>
                <div class="cfg-select">
                    <select required>
                        <option value="" disabled selected hidden>Выберите</option>
                        <option value="1">Выберите размерность пространства</option>
                        <option value="2">Выберите размерность пространства</option>
                    </select>
                </div>
            </div>
            <div class="cfg-select-inf-item">
                <div class="cfg-lbl">
                    Границы пространства
                    <div class="cfg-hint-bx">
                        <div class="cfg-opener"><img
                                src="<?php echo $server_url; ?>/construct/static/images/info-circle.svg" alt="">
                        </div>
                        <div class="cfg-frame">Всплывающая информация по данному полю</div>
                    </div>
                </div>
                <input type="text" placeholder="Число">
                <div class="cfg-sel-over">
                    <div class="cfg-lbl">
                        Ед. изм.
                        <div class="cfg-hint-bx">
                            <div class="cfg-opener"><img
                                    src="<?php echo $server_url; ?>/construct/static/images/info-circle.svg" alt="">
                            </div>
                            <div class="cfg-frame">Всплывающая информация по данному полю</div>
                        </div>
                    </div>
                    <div class="cfg-select">
                        <select required>
                            <option value="" disabled selected hidden>Выберите</option>
                            <option value="1">Ед. Изм.</option>
                        </select>
                    </div>
                </div>
            </div>
            <div class="cfg-ch-group">
                <label class="cfg-check-item">
                    <input type="checkbox">
                    <span>
                        Учитывать гравитационные взаимодействия
                        <div class="cfg-hint-bx">
                            <div class="cfg-opener"><img
                                    src="<?php echo $server_url; ?>/construct/static/images/info-circle.svg" alt="">
                            </div>
                            <div class="cfg-frame">Всплывающая информация по данному полю</div>
                        </div>
                    </span>
                </label>
                <label class="cfg-check-item">
                    <input type="checkbox">
                    <span>
                        Учитывать мой стаж работы
                        <div class="cfg-hint-bx">
                            <div class="cfg-opener"><img
                                    src="<?php echo $server_url; ?>/construct/static/images/info-circle.svg" alt="">
                            </div>
                            <div class="cfg-frame">Всплывающая информация по данному полю</div>
                        </div>
                    </span>
                </label>
            </div>
            <div class="cfg-select-inf-item">
                <div class="cfg-lbl">
                    Коэффицент восстановления
                    <div class="cfg-hint-bx">
                        <div class="cfg-opener"><img
                                src="<?php echo $server_url; ?>/construct/static/images/info-circle.svg" alt="">
                        </div>
                        <div class="cfg-frame">Всплывающая информация по данному полю</div>
                    </div>
                </div>
                <input type="text" placeholder="Число">
            </div>
            <div class="cfg-back-bx">
                <a href="#cfg-page-1" class="cfg-btn cfg-btn-framed">Назад</a>
                <a href="#cfg-page-3" class="cfg-btn cfg-btn-next">Далее</a>
            </div>
        </div>
    </div>
</div>


<div class="cfg-container" id="cfg-page-3">

    <div class="cfg-table-bx">
        <div class="cfg-h3">Создание объектов задачи</div>

        <div class="cfg-step-bx">

            <div class="cfg-table-wrap">
                <table class="cfg-table cfg-table-red">
                    <tr>
                        <td>Название</td>
                        <td>Тип</td>
                        <td>Состояние</td>
                        <td>Цвет</td>
                        <td></td>
                    </tr>
                    <tr class="cfg-filled">
                        <td>Солнце</td>
                        <td>Звезда</td>
                        <td>Горячее</td>
                        <td>
                            <div class="cfg-color cfg-color1"></div>
                        </td>
                        <td>
                            <div class="cfg-red cfg-pop-form-opener"><img
                                    src="<?php echo $server_url; ?>/construct/static/images/red-bt.svg" alt=""></div>
                        </td>
                    </tr>
                    <tr class="cfg-filled">
                        <td>Солнце</td>
                        <td>Звезда</td>
                        <td>Горячее</td>
                        <td>
                            <div class="cfg-color cfg-color2"></div>
                        </td>
                        <td>
                            <div class="cfg-red cfg-pop-form-opener"><img
                                    src="<?php echo $server_url; ?>/construct/static/images/red-bt.svg" alt=""></div>
                        </td>
                    </tr>
                    <tr class="cfg-filled">
                        <td>Солнце</td>
                        <td>Звезда</td>
                        <td>Горячее</td>
                        <td>
                            <div class="cfg-color cfg-color3"></div>
                        </td>
                        <td>
                            <div class="cfg-red cfg-pop-form-opener"><img
                                    src="<?php echo $server_url; ?>/construct/static/images/red-bt.svg" alt=""></div>
                        </td>
                    </tr>
                </table>
                <a href="" class="cfg-plus-bt"><img src="<?php echo $server_url; ?>/construct/static/images/plus2.svg"
                        alt=""></a>
            </div>

            <div class="cfg-back-bx">
                <a href="#cfg-page-2" class="cfg-btn cfg-btn-framed">Назад</a>
                <button class="cfg-btn cfg-btn-next">Сохранить</button>
            </div>
        </div>

    </div>

</div>


<div class="cfg-overlay"></div>

<div class="cfg-popup cfg-popup-form">
    <div class="cfg-close"><img src="<?php echo $server_url; ?>/construct/static/images/close.svg" alt=""></div>

    <form action="" class="cfg-pop-form">
        <div class="cfg-h2">Настройки параметров объекта</div>
        <div class="cfg-wr cfg-wr1">
            <div class="cfg-select-inf-item cfg-bx1">
                <div class="cfg-lbl">
                    Название объекта
                    <div class="cfg-hint-bx">
                        <div class="cfg-opener"><img
                                src="<?php echo $server_url; ?>/construct/static/images/info-circle.svg" alt=""></div>
                        <div class="cfg-frame">Всплывающая информация по данному полю</div>
                    </div>
                </div>
                <input type="text" placeholder="Текст">
            </div>
            <div class="cfg-select-inf-item cfg-bx2 cfg-select-inf-item-color">
                <div class="cfg-lbl">
                    Цвет объекта
                    <div class="cfg-hint-bx">
                        <div class="cfg-opener"><img
                                src="<?php echo $server_url; ?>/construct/static/images/info-circle.svg" alt=""></div>
                        <div class="cfg-frame">Всплывающая информация по данному полю</div>
                    </div>
                </div>
                <div class="color"></div>
                <input type="text" value="#60C64D">
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
                <select>
                    <option value="">Выберите</option>
                    <option value="">Единичный объект</option>
                </select>
            </div>
        </div>
        <div class="cfg-select-inf-item">
            <div class="cfg-lbl">
                Координаты центра
                <div class="cfg-hint-bx">
                    <div class="opener"><img src="<?php echo $server_url; ?>/construct/static/images/info-circle.svg"
                            alt=""></div>
                    <div class="frame">Всплывающая информация по данному полю</div>
                </div>
            </div>
            <div class="cfg-xyz-bx">
                <div class="cfg-bx">
                    <div class="cfg-let">X</div>
                    <input type="text" placeholder="Число">
                </div>
                <div class="cfg-bx">
                    <div class="cfg-let">Y</div>
                    <input type="text" placeholder="Число">
                </div>
                <div class="cfg-bx">
                    <div class="cfg-let">Z</div>
                    <input type="text" placeholder="Число">
                </div>
            </div>
        </div>

        <div class="cfg-select-inf-item">
            <div class="cfg-line-bx">
                <div class="cfg-lbl">
                    Ед. измерения
                    <div class="cfg-hint-bx">
                        <div class="cfg-opener"><img
                                src="<?php echo $server_url; ?>/construct/static/images/info-circle.svg" alt=""></div>
                        <div class="cfg-frame">Всплывающая информация по данному полю</div>
                    </div>
                </div>
                <div class="cfg-select">
                    <select>
                        <option value="">Выберите</option>
                        <option value="">Единичный объект</option>
                    </select>
                </div>
            </div>
        </div>


        <div class="cfg-select-inf-item">
            <div class="cfg-lbl">
                Заряд
                <div class="cfg-hint-bx">
                    <div class="cfg-opener"><img
                            src="<?php echo $server_url; ?>/construct/static/images/info-circle.svg" alt=""></div>
                    <div class="cfg-frame">Всплывающая информация по данному полю</div>
                </div>
            </div>
            <input type="text" placeholder="Число">
            <div class="cfg-sel-over">
                <div class="cfg-lbl">
                    Ед. изм.
                    <div class="cfg-hint-bx">
                        <div class="cfg-opener"><img
                                src="<?php echo $server_url; ?>/construct/static/images/info-circle.svg" alt=""></div>
                        <div class="cfg-frame">Всплывающая информация по данному полю</div>
                    </div>
                </div>
                <div class="cfg-select">
                    <select>
                        <option value="">Выберите</option>
                        <option value="">Ед. Изм.</option>
                    </select>
                </div>
            </div>
        </div>

        <div class="cfg-select-inf-item">
            <div class="cfg-lbl">
                Масса
                <div class="cfg-hint-bx">
                    <div class="cfg-opener"><img
                            src="<?php echo $server_url; ?>/construct/static/images/info-circle.svg" alt=""></div>
                    <div class="cfg-frame">Всплывающая информация по данному полю</div>
                </div>
            </div>
            <input type="text" placeholder="Число">
            <div class="cfg-sel-over">
                <div class="cfg-lbl">
                    Ед. изм.
                    <div class="cfg-hint-bx">
                        <div class="cfg-opener"><img
                                src="<?php echo $server_url; ?>/construct/static/images/info-circle.svg" alt=""></div>
                        <div class="cfg-frame">Всплывающая информация по данному полю</div>
                    </div>
                </div>
                <div class="cfg-select">
                    <select>
                        <option value="">Выберите</option>
                        <option value="">Ед. Изм.</option>
                    </select>
                </div>
            </div>
        </div>

        <div class="cfg-select-inf-item">
            <div class="cfg-lbl">
                Радиус
                <div class="cfg-hint-bx">
                    <div class="cfg-opener"><img
                            src="<?php echo $server_url; ?>/construct/static/images/info-circle.svg" alt=""></div>
                    <div class="cfg-frame">Всплывающая информация по данному полю</div>
                </div>
            </div>
            <input type="text" placeholder="Число">
            <div class="cfg-sel-over">
                <div class="cfg-lbl">
                    Ед. изм.
                    <div class="cfg-hint-bx">
                        <div class="cfg-opener"><img
                                src="<?php echo $server_url; ?>/construct/static/images/info-circle.svg" alt=""></div>
                        <div class="cfg-frame">Всплывающая информация по данному полю</div>
                    </div>
                </div>
                <div class="cfg-select">
                    <select>
                        <option value="">Выберите</option>
                        <option value="">Ед. Изм.</option>
                    </select>
                </div>
            </div>
        </div>

        <div class="cfg-bt-wr">
            <a href="" class="cfg-btn cfg-btn-framed">Удалить</a>
            <input type="submit" class="cfg-btn" value="Сохранить">
        </div>
    </form>
</div>