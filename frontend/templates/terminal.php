
<link rel="preconnect" href="https://fonts.googleapis.com">
<link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
<link href="https://fonts.googleapis.com/css2?family=Open+Sans:wght@400;700&display=swap" rel="stylesheet">
<link rel="stylesheet" href="<?php echo $server_url; ?>/construct/static/css/config.css">
<script type="text/javascript" src="https://cdnjs.cloudflare.com/ajax/libs/jquery/3.5.0/jquery.min.js"></script>
<script type="text/javascript" src="https://cdnjs.cloudflare.com/ajax/libs/FileSaver.js/2.0.0/FileSaver.min.js"></script>
<script type="text/javascript" src="<?php echo $server_url; ?>/construct/static/js/terminal.js"></script>

<div class="cfg-container" style="display:block;opacity:0;max-width:none;">

	<div class="cfg-table-bx">
		<div style="display:flex;align-items:center;margin-bottom:20px">
			<div class="cfg-h3" style="margin:0;">Создание объектов задачи</div>
			<input type="file" id="cfg-file-input" style="display:none;" accept=".json" />
            <img class="cfg-plus-bt" onclick="$('#cfg-file-input').trigger('click');" style="margin:0;margin-left:auto;"src="<?php echo $server_url; ?>/construct/static/images/plus2.svg">
		</div>
		<div class="cfg-table-wrap">
			<table class="cfg-table" id="cfg-table">
				<tr>
					<td>№</td>
					<td>Сценарий</td>
					<td>Тип проблемы</td>
					<td>Состояние</td>
					<td>Создана</td>
					<td>Завершена</td>
					<td></td>
				</tr>
			</table>
		</div>
	</div>

</div>

<div class="cfg-overlay"></div>
<div class="cfg-popup cfg-popup-video" id="cfg-result"></div>
