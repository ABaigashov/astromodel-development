from traceback import format_exc as _traceback
from configurator import Configurator
import json, os, pickle, subprocess
from threading import Thread
from blessed import Terminal


IS_SERVER = os.getcwd() == '/usr/src'

class Dispatcher:

	def __init__(self, config, output, job=None):
		with open(config, 'rb') as f:
			parameters = json.load(f)
		problem, configuration, animation_length = Configurator.parse_parameters(parameters)
		self.animation_length = animation_length
		self.configuration = configuration
		self.problem = problem
		self.output = output
		self.job = job

	def init(self):
		self.process = subprocess.Popen("/bin/bash",
			stderr=subprocess.STDOUT,
			stdout=subprocess.PIPE,
			stdin=subprocess.PIPE,
			env=os.environ.copy(),
			encoding='utf-8',
			shell=True
		)
		if IS_SERVER:
			self.process.stdin.write('export WORKON_HOME=$pwd/enviroments\n')
			self.process.stdin.write('source /usr/local/bin/virtualenvwrapper.sh\n')
			self.process.stdin.write(f'workon {self.problem}\n')
			self.process.stdin.write(f'python3 -u -B modeling_module/physical_problems/{self.problem}/main.py {self.configuration} {self.output}\n')
		else:
			self.process.stdin.write(f'python3 -u -B main.py {self.configuration} {self.output}\n')
		self.process.stdin.write('exit\n')
		self.process.stdin.flush()

	def run(self):
		output = ''
		counter = 0
		if not IS_SERVER:
			term = Terminal()
		while True:
			try:
				symbol = self.process.stdout.read(1)
				if not symbol and self.process.poll() is not None:
					break
				if IS_SERVER:
					if ']' in output[-1:] and symbol == '\n':
						self.job.progress = min(0.99, counter / self.animation_length)
						counter += 1
				else:
					if ']' in output[-1:] and symbol == '\n':
						print(end=symbol + term.move_up(1))
					elif output and '\n' in output[-1:] and symbol != 't':
						print(end=term.move_down(1) + symbol)
					else:
						print(end=symbol)
				output += symbol
				self.process.stdout.flush()
			except KeyboardInterrupt:
				exit()
			except Exception as e:
				return repr(e)
		return output.split('\n')[-2]

class JobExecutor(Thread):

	def __init__(self, name):
		super().__init__()
		self.name = name

	def Init(self, job, config, output):
		self.job = job
		self.output = output
		self.config = config
		self.path_to_result = os.path.join(output, 'error.txt')
		try:
			self.dispatcher = Dispatcher(config, output, job)
			self.dispatcher.init()
		except:
			print(_traceback())
			with open(self.path_to_result, 'wb') as logfile:
				logfile.write(_traceback().encode())
			self.job.progress = -1

	def run(self):
		try:
			temp = self.dispatcher.run()
			if not os.path.isfile(temp):
				raise eval(temp)
			self.path_to_result = temp
			self.job.progress = 1
		except:
			print(_traceback())
			with open(self.path_to_result, 'wb') as logfile:
				logfile.write(_traceback().encode())
			self.job.progress = -1
