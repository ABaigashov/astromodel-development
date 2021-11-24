
# This is main problem file. It has special name
# '__init__.py' and you SHOULDN'T RENAME it. The file
# contains main solver object.


# Import LOCAL python files with 'GlobalInteraction' object
# and specifiend models representation of problem solution.
from solver import Cosmology_data, Cosmology_calculus, Visualization, Model_Var
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
		self.config = config
		self.model1 = config.cosmological_components
		self.models = []
		self.model_SNE = Cosmology_data(config.name_SNE, [int(config.row_1),int(config.row_2),int(config.row_3)])
		self.model_H = Cosmology_data(config.name_Hubble,[int(config.row_4),int(config.row_5),int(config.row_6)])
		self.model_SNE.Data_loader()
		self.model_H.Data_loader()

		if config.task =="0":
			# saving current problem model with incomming parameters
			for model in config.cosmological_components:
				self.models.append(Cosmology_calculus(config, model, self.model_SNE, self.model_H))

			for model in self.models:
				model.mu_diagram()
				model.integration()
				model.hubble_versus_z()
				model.mu_versus_z()

			self.GRAPH = Visualization(self.models, output, job)

		else:
			for model in config.cosmological_components:
				print(model)
				omega_m = float(model.omega_m)+0.1
				omega_d = float(model.omega_d)-0.1
				omega_r = float(model.omega_r)
				equation_d = model.equation_d
				title = '0'
				model_var = Model_Var(str(omega_d), str(omega_m), str(omega_r), equation_d, title)

				print(omega_m,omega_d,omega_r,equation_d)
				self.models.append(Cosmology_calculus(config, model_var, self.model_SNE, self.model_H))

			for model in self.models:
				model.mu_diagram()
				model.integration()
				model.hubble_versus_z()
				model.mu_versus_z()

			self.GRAPH = Visualization(self.models, output, job)








	# run method
	# must ALWAYS return path to rendered file
	def run(self):

		# render file and return path
		return self.GRAPH.graphics(self.config)
