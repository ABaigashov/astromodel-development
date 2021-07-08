from traceback import format_exc as _traceback, format_exception as error_text
import json, os, pickle, subprocess, sys
from configurator import Configurator
from threading import Thread
from uuid import uuid4

IS_WINDOWS = os.name == 'nt'
python = 'python' if IS_WINDOWS else 'python3'
command = 'cmd' if IS_WINDOWS else '/bin/bash'
IS_SERVER = os.getcwd() == '/usr/src'

class Dispatcher:

	def __init__(self, config, output, job=None):
		with open(config, 'rb') as f:
			parameters = json.load(f)
		problem, configuration = Configurator.parse_parameters(parameters)
		self.configuration = configuration
		self.problem = problem
		self.output = output
		self.job = job
		self.load_backcall()

	def load_backcall(self):
		if IS_SERVER:
			sys.path.append(os.path.join('.', 'modeling_module', 'physical_problems', self.problem))
		else:
			sys.path.append(os.path.join(os.path.dirname(__file__), 'modeling_module', 'physical_problems', self.problem))
		try:
			from backcall import value
			self.backcall_count = value(Configurator(self.configuration))
		except Exception as e:
			self.backcall_count = 1
		sys.path.pop()

	def init(self):
		self.process = subprocess.Popen(command,
			stderr=subprocess.STDOUT,
			stdout=subprocess.PIPE,
			stdin=subprocess.PIPE,
			env=os.environ.copy(),
			shell=True
		)
		self.exitkey = '[exitkey:' + uuid4().hex + ']'
		if IS_SERVER:
			venv = os.path.join('.', 'enviroments', self.problem, 'bin', 'python3')
			self.command = f'{venv if os.path.exists(venv) else python} \
				-u -B modeling_module/physical_problems/{self.problem}/main.py \
				{pickle.dumps(self.configuration).hex()} {self.output} {self.exitkey}\n'
		else:
			self.command = f'{python} -u -B main.py {pickle.dumps(self.configuration).hex()} {self.output} {self.exitkey}\n'
		self.process.stdin.write(self.command.encode())
		self.process.stdin.write(b'exit\n')
		self.process.stdin.flush()

	def run(self):
		output = ''
		counter = 0
		SYMBOLS_STACK = b''
		flag = end = False
		try:
			while True:
				symbol = self.process.stdout.read(1)
				try:
					if symbol == b'\b':
						if IS_SERVER:
							self.job.progress = min(0.9999, counter / abs(self.backcall_count))
							counter += 1
						continue
					SYMBOLS_STACK += symbol
					symbol = SYMBOLS_STACK.decode()
					SYMBOLS_STACK = b''
				except:
					if IS_WINDOWS:
						output = ''
						SYMBOLS_STACK = b''
					continue
				if not symbol and self.process.poll() is not None:
					break
				if IS_WINDOWS and self.command and symbol == self.command[0]:
					flag = True
					self.command = self.command[1:]
				if not IS_WINDOWS or (flag and not len(self.command)):
					if not IS_SERVER:
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
