import os, sys, time
import json, codecs, shutil
import asyncio
import hashlib
import websockets
import traceback
import random
from threading import Thread
from database.Job import Job
from database.JobNode import JobNode

sys.path.append(os.path.join(os.getcwd(), 'modeling_module'))

from jobexecutor import JobExecutor

def get_path(job):
	path = "/var/astromodel/current_" + str(job.uid) + "/";
	if not os.path.exists(path):
		os.mkdir(path)
	return path

class Request:
	def __init__(self, command, **kwargs):
		self.command = command
		self.args = kwargs or None

	def toJSON(self):
		return json.dumps(self, cls=RequestEncoder)


class RequestEncoder(json.JSONEncoder):
	def default(self, obj):
		if isinstance(obj, Request):
			return obj.__dict__
		return json.JSONEncoder.default(self, obj)


class AstroNode:

	def __init__(self, host, name):
		self.host = host
		self.name = name
		self.job = None
		self.jobnode = None
		self.connected = False
		self.executor_thread = None

	def update_state(self):
		if self.job is None:
			self.jobnode.state = 'idle'
			self.jobnode.job_uid = 0
			self.jobnode.progress = 0
			self.jobnode.problem = os.environ['PROBLEM']
			# self.jobnode.extention = ''
		else:
			self.jobnode.job_uid = self.job.uid
			self.jobnode.progress = self.job.progress

	async def wsstatus(self):
		if self.jobnode is None:
			request = Request("handshake", name=self.name)
			return request.toJSON()
		else:
			self.update_state()
			print('Update node status')
			request = Request("nodestatus", **self.jobnode.__dict__)
			return request.toJSON()

	async def wsconsumer(self, websocket):
		try:
			message = await asyncio.wait_for(websocket.recv(), 1)
		except asyncio.TimeoutError:
			message = None
		except websockets.exceptions.ConnectionClosedError as e:
			raise e
		return message

	async def wsstart(self):
		#try:
		self.executor_thread = None
		async with websockets.connect(self.host, max_size=1048576*50) as websocket:
			print('Connection established')
			self.connected = True
			while self.connected:
				if self.executor_thread is None:
					if self.jobnode is not None and self.job is not None:
						self.jobnode.state = 'execution'
						path_result = get_path(self.job)
						path_config = path_result + "config.json"
						self.executor_thread = JobExecutor("JobExecutor")
						self.executor_thread.Init(self.job, path_config, path_result)
						self.executor_thread.start()
				elif self.job.progress == 1:
					self.executor_thread.join()
					self.path_to_result = self.executor_thread.path_to_result
					self.executor_thread = None
					self.jobnode.extention = '.' + self.path_to_result.split('.')[-1]
					message = await self.wsstatus()
					await websocket.send(message)

					self.jobnode.state = 'complete'
					file_result = open(self.path_to_result, "rb")
					data_result = file_result.read()
					file_result.close()

					await websocket.send(data_result)
					print('Job done')

				elif self.job.progress == -1:
					self.executor_thread.join()
					self.path_to_result = self.executor_thread.path_to_result
					
					self.executor_thread = None
					self.jobnode.extention = '.txt'

					message = await self.wsstatus()
					await websocket.send(message)

					self.jobnode.state = 'failure'
					file_result = open(self.path_to_result, "rb")

					data_result = file_result.read()

					file_result.close()
					await websocket.send(data_result)
					print('Job aborted')

				message = await self.wsstatus()
				await websocket.send(message)
				await asyncio.sleep(1)

				message = await self.wsconsumer(websocket)
				if message is not None:
					response = json.loads(message)
					if "command" in response.keys():
						if response["command"] == "register":
							self.jobnode = JobNode()
							self.jobnode.FromData(**response["args"])
							self.jobnode.problem = os.environ['PROBLEM']
							print("JobNode has been registered")
						elif response["command"] == "execute":
							if 'job' in response['args']:
								self.job = Job()
								self.job.FromData(**response["args"]['job'])

							path_result = get_path(self.job)
							path_config = path_result + "config.json"

							if 'config' in response['args']:
								f = open(path_config,"w+")
								f.write(json.dumps(response["args"]['config']))
								f.close()
							else:
								if os.path.exists(path_config):
									os.remove(path_config)
							print("Job has been accepted")
						elif response["command"] == "release":
							path_result = get_path(self.job)
							if os.path.exists(path_result):
								shutil.rmtree(path_result)

							job_task = None
							self.job = None
							self.update_state()
							print("JobNode has been released")
		#except Exception as ex:
		#    print(str(ex))
	#        print('Connection failed')

def random_uniq(length):
	return random_string(length) + '-' + random_string(length)

def random_string(length):
	alphabet = 'abcdefghijklmnopqrstuvwxyz'
	return ''.join((random.choice(alphabet) for i in range(length)))

try:
	client = AstroNode('ws://wsserver:8080', random_uniq(3))
	loop = asyncio.get_event_loop()
	loop.run_until_complete(client.wsstart())
	loop.run_forever()
except KeyboardInterrupt:
	print("Client is shuttong down...")
