
from contextlib import redirect_stdout
from configurator import Configurator
from traceback import format_exc
from threading import Thread
from uuid import uuid4
import json, os, sys


class JobExecutor(Thread):

	def __init__(self, name):
		super().__init__()
		self.problem = os.environ['PROBLEM']
		sys.path.append(os.path.join('.', 'modeling_module', 'physical_problems', self.problem))
		exec(f'from physical_problems.{self.problem} import Model', {}, self.__dict__)
		self.name = name

	def Init(self, job, configuration, output):
		self.job = job
		self.path_to_result = os.path.join(output, 'error.txt')
		try:
			config = Configurator(configuration)
			self.model = self.Model()
			self.model.init(config, os.path.join(output, uuid4().hex[:16]), job)
		except:
			print(format_exc())
			with open(self.path_to_result, 'wb') as logfile:
				logfile.write(format_exc().encode())
			self.job.progress = -1

	def run(self):
		try:
			with open(os.devnull, 'w') as devnull:
				with redirect_stdout(devnull):
					self.path_to_result = self.model.run()
				self.job.progress = 1
		except:
			print(format_exc())
			with open(self.path_to_result, 'wb') as logfile:
				logfile.write(format_exc().encode())
			self.job.progress = -1
