
# This is main problem file. It has special name
# '__init__.py' and you SHOULDN'T RENAME it. The file
# contains main solver object.


# Import LOCAL python files with 'GlobalInteraction' object
# and specifiend models representation of problem solution.
from solver import Task_maker, Cosmology_data, Cosmology_calculus
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


		Task = Task_maker(config)
		# saving current problem model with incomming parameters

		self.model_SNE = Cosmology_data(Task.name_SNE, Task.row_SNE)
		self.model_H = Cosmology_data(Task.name_Hubble, Task.row_Hubble)
		self.model_SNE.Data_loader()
		self.model_H.Data_loader()
		self.calculations = Cosmology_calculus(Task, self.model_SNE, self.model_H)

	# run method
	# must ALWAYS return path to rendered file
	def run(self):
		self.calculations.mu_diagram()
		self.calculations.integration()
		self.calculations.visualization()

		# render file and return path
		return 1
