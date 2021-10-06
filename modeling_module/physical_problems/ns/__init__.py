
# This is main problem file. It has special name
# '__init__.py' and you SHOULDN'T RENAME it. The file
# contains main solver object.


# Import LOCAL python files with 'GlobalInteraction' object
# and specifiend models representation of problem solution.
from graphical_output import Model_of_representation
from task_maker import Task_maker
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
		Task.EoS_append()
		Task.Task_append()

		if config.calculate_mass_profile == True:
			self.task2 = "yes"
		else:
			self.task2 = "no"

		if config.calculate_mass_in_interval == True:
			self.task1 = "yes"
		else:
			self.task1 = "no"

		self.model = Model_of_representation(Task.task_descriptor, Task.EoS_descriptor)

	# run method
	# must ALWAYS return path to rendered file
	def run(self):
		if self.task2=="yes":
			self.model.work()
		if self.task1=="yes":
			self.model.work_2()
		print("Вычисления завершены, посмотрите результаты в папке results")
		# render file and return path
		return 0
