
# This is main problem file. It has special name
# '__init__.py' and you SHOULDN'T RENAME it. The file
# contains main solver object.


# Import LOCAL python files with 'GlobalInteraction' object
# and specifiend models representation of problem solution.
from solver import Task_maker, BVP_solver
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
		if self.config.mass_object:
			self.potential = PointsPotential(self.config, output, job)
		else:
			Task = Task_maker(self.config)
			# saving current problem model with incomming parameters

			self.model = BVP_solver(Task, output)

	# run method
	# must ALWAYS return path to rendered file
	def run(self):
		if self.config.mass_object:
			return self.potential.start()
		else:
			# print("Вычисления завершены, посмотрите результаты в папке results")
			# render file and return path
			return self.model.Solving_eq()
