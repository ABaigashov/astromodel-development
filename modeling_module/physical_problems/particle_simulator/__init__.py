
from interaction_creator import GlobalInteraction
from models import _model_from_config

class Model:

	def init(self, config, output, job):
		astro_object = GlobalInteraction(config)
		self.model = _model_from_config(config, astro_object, output, job)

	def run(self):
		return self.model.render()
