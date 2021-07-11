from traceback import format_exc as _traceback, format_exception as error_text
import json, os, pickle, subprocess, sys
from configurator import Configurator
from threading import Thread
from uuid import uuid4


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
		sys.path.append(os.path.join('.', 'modeling_module', 'physical_problems', self.problem))
		try:
			from backcall import value
			self.backcall_count = value(Configurator(self.configuration))
		except Exception as e:
			self.backcall_count = 1
		sys.path.pop()

	def init(self):
		self.process = subprocess.Popen('/bin/bash',
			stderr=subprocess.STDOUT,
			stdout=subprocess.PIPE,
			stdin=subprocess.PIPE,
			env=os.environ.copy(),
			shell=True
		)
		self.exitkey = '[exitkey:' + uuid4().hex + ']'
		venv = os.path.join('.', 'enviroments', self.problem, 'bin', '')
		self.command = f'{venv if os.path.exists(venv) else ''}python3 -u -B \
			modeling_module/physical_problems/{self.problem}/main.py \
			{pickle.dumps(self.configuration).hex()} {self.output} {self.exitkey}\n'
		self.process.stdin.write(self.command.encode())
		self.process.stdin.write(b'exit\n')
		self.process.stdin.flush()

	def run(self):
		output = ''
		counter = 0
		SYMBOLS_STACK = b''
		end = False
		try:
			while True:
				symbol = self.process.stdout.read(1)
				if symbol == b'\b':
					self.job.progress = min(0.9999, counter / self.backcall_count)
					counter += 1
					continue
				try:
					SYMBOLS_STACK += symbol
					symbol = SYMBOLS_STACK.decode()
					SYMBOLS_STACK = b''
				except:
					continue
				if not symbol and self.process.poll() is not None:
					break
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
		try:
			return output.split('\n')[-2].split(self.exitkey)[1][1:]
		except:
			return

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
			self.path_to_result = result
			if isinstance(result, tuple):
				error = ''.join(error_text(*result))
				print(error)
				with open(self.path_to_result, 'wb') as logfile:
					logfile.write(error.encode())
				self.job.progress = -1
				return
			self.job.progress = 1
		except:
			print(_traceback())
			with open(self.path_to_result, 'wb') as logfile:
				logfile.write(_traceback().encode())
			self.job.progress = -1
