

# This is main problem file. It has special name
# '__init__.py' and you SHOULDN'T RENAME it. The file
# contains main solver object.


# Import LOCAL python files with 'GlobalInteraction' object
# and specifiend models representation of problem solution.
from interaction_creator import GlobalInteraction
from models import _model_from_config


# This is the main class. It has special name and
# in your '__init__.py' file it MUST have same name
# and same methods.
class Model:

	# init method
	# Arguments :self: class object
	#         :config: problem configuration (instance of 'Configurator')
	#         :output: string with output path (without extention)
	#            :job: job instance (uses for 'job.process' monipulations)
	def init(self, config, output, job):

		# creating object of 'GlobalInteraction' class
		astro_object = GlobalInteraction(config)

		# saving current problem model with incomming parameters
		self.model = _model_from_config(config, astro_object, output, job)

	# run method
	# must ALWAYS return path to rendered file
	def run(self):

		# render file and return path
		return self.model.render()
