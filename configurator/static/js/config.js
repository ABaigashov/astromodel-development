jQuery(document).ready(function($){
	
	$('.pop-video').on('click', function(event){
		event.preventDefault()
		$('.overlay, .popup-video').fadeIn();
	})
	$('.pop-form-opener').on('click', function(event){
		event.preventDefault()
		$('.overlay, .popup-form').fadeIn();
	})
	$('.overlay, .popup .close').on('click', function(event){
		$('.overlay, .popup').fadeOut();
	})

})