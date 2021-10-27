
# This is main problem file. It has special name
# '__init__.py' and you SHOULDN'T RENAME it. The file
# contains main solver object.


# Import LOCAL python files with 'GlobalInteraction' object
# and specifiend models representation of problem solution.
from solver import Task_maker, Cosmology_data, Cosmology_calculus, Visualization
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


		self.Task = Task_maker(config)
		# saving current problem model with incomming parameters
		self.models = []
		self.models_1 = []
		self.model_SNE = Cosmology_data(self.Task.name_SNE, self.Task.row_SNE)
		self.model_H = Cosmology_data(self.Task.name_Hubble, self.Task.row_Hubble)
		self.model_SNE.Data_loader()
		self.model_H.Data_loader()
		for model in self.Task.cosmological_components:
			self.models.append(Cosmology_calculus(config, model, self.model_SNE, self.model_H))

		for model in self.models:
			model.mu_diagram()
			model.integration()
			model.hubble_versus_z()
			model.chi_square_hubble()
			# print(model.chi_square_H, model.H_opt)
			self.models_1.append(model)

		self.GRAPH = Visualization(self.models_1, output)

	# run method
	# must ALWAYS return path to rendered file
	def run(self):

		# render file and return path
		return self.GRAPH.graphics(self.Task)
