from traceback import format_exc as _traceback, format_exception as error_text
import json, os, pickle, subprocess, sys
from configurator import Configurator
from threading import Thread
from blessed import Terminal
from uuid import uuid4

python = 'python' if os.name == 'nt' else 'python3'
command = 'cmd' if os.name == 'nt' else '/bin/bash'
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
		self.process = subprocess.Popen(command,
			stderr=subprocess.STDOUT,
			stdout=subprocess.PIPE,
			stdin=subprocess.PIPE,
			env=os.environ.copy(),
			encoding='utf-8',
			shell=True
		)
		self.exitkey = '[exitkey:' + uuid4().hex + ']'
		if IS_SERVER:
			self.process.stdin.write('export WORKON_HOME=$pwd/enviroments\n')
			self.process.stdin.write('source /usr/local/bin/virtualenvwrapper.sh\n')
			self.process.stdin.write(f'workon {self.problem}\n')
			self.command = f'{python} -u -B modeling_module/physical_problems/{self.problem}/main.py {self.configuration} {self.output} {self.exitkey}\n'
		else:
			self.command = f'{python} -u -B main.py {self.configuration} {self.output} {self.exitkey}\n'
		self.process.stdin.write(self.command)
		self.process.stdin.write('exit\n')
		self.process.stdin.flush()

	def run(self):
		output = ''
		counter = 0
		flag = end = False
		if not IS_SERVER:
			term = Terminal()
		while True:
			try:
				try:
					symbol = self.process.stdout.read(1)
				except:
					output = ''
				if self.command and symbol == self.command[0]:
					flag = True
					self.command = self.command[1:]
				if not symbol and self.process.poll() is not None:
					break
				if flag and not len(self.command):
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
					if output.endswith(self.exitkey):
						end = True
					if end and symbol == '\n':
						break
			except KeyboardInterrupt:
				exit()
			except Exception:
				return sys.exc_info()
		return output.split('\n')[-2].split(self.exitkey)[1][1:]

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
			result = self.dispatcher.run()
			if isinstance(result, tuple):
				error = ''.join(error_text(*result))
				print(error)
				with open(self.path_to_result, 'wb') as logfile:
					logfile.write(error.encode())
				self.job.progress = -1
				return
			self.path_to_result = result
			self.job.progress = 1
		except:
			print(_traceback())
			with open(self.path_to_result, 'wb') as logfile:
				logfile.write(_traceback().encode())
			self.job.progress = -1
