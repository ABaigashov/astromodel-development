import os, sys, uuid
import asyncio, aiopg, websockets, json, codecs
from database.Database import Database
from database.JobNode import JobNode
from database.Job import Job

path_root = "/var/astromodel/"


class JobNodeClient:
	def __init__(self, websocket):
		self.websocket = websocket

	async def set_node(self, name):
		self.jobnode = await Database.JobNodeLoadByName(name)
		if self.jobnode is None:
			self.jobnode = await Database.JobNodeCreate(str(uuid.uuid4()), name, 0)

		self.jobnode.state = 'none'
		self.jobnode.extention = ''
		self.jobnode.job_uid = 0
		self.jobnode.progress = 0

	async def update(self, **kwargs):
		self.jobnode.__dict__.update(kwargs)
		await self.jobnode.Update()


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


class AstroServer:

	def __init__(self):
		self.wsclients = set()
		self.job_busy = set()

	async def find(self, websocket):
		node = None
		for c in self.wsclients:
			if c.websocket == websocket:
				node = c
				break
		return node

	async def register(self, node):
		self.wsclients.add(node)

	async def unregister(self, websocket):
		node = await self.find(websocket)
		if node is not None:
			self.wsclients.remove(node)

	async def wsconsumer(self, node):
		try:
			message = await asyncio.wait_for(node.websocket.recv(), 1)
			if not isinstance(message, str):
				job = await Database.JobLoad(node.jobnode.job_uid)
				job_result = job.path + 'result' + node.jobnode.extention
				print("Result path: " + str(job_result))
				result_file = open(job_result, "wb")
				result_file.write(message)
				result_file.close()
				message = None
		except asyncio.TimeoutError:
			message = None
		except websockets.exceptions.ConnectionClosedError as e:
			raise e
		return message

	async def readconfig(self, path):
		config = None
		if os.path.exists(path):
			file = codecs.open(path, "r", "utf-8")
			config = json.loads(file.read())
			file.close()
		return config

	async def connect(self, websocket, path):
		print('Connection established:', path)
		if('/health' == path):
			return

		try:
			node = JobNodeClient(websocket)
			while True:
				request = None
				if node is not None:
					if 'jobnode' in node.__dict__ and node.jobnode is not None:
						if node.jobnode.state == 'complete':
							if node.jobnode.job_uid > 0:
								job = await Database.JobLoad(node.jobnode.job_uid)
								job.extention = node.jobnode.extention
								await job.Finish()
								print("Job in done: " + str(node.jobnode.job_uid))

							request = {}
							request = Request("release", **request)
							print("Job result has been retrieved")
							await websocket.send(request.toJSON())
						elif node.jobnode.state == 'failure':
							job = await Database.JobLoad(node.jobnode.job_uid)
							job.extention = node.jobnode.extention
							await job.Error()
							await Database.JobUpdate(job)

							job_result = job.path + 'result' + job.extention
							# if os.path.exists(job_result):
							# 	os.rename(job_result, job.path + "error.txt")

							print("Job error: " + str(node.jobnode.progress))
							request = {}
							request = Request("release", **request)
							print("Job result has been retrieved")
							await websocket.send(request.toJSON())
						elif node.jobnode.state == 'idle':
							job_awaits = await Database.JobLoadAwaiting()
							if job_awaits is not None:
								config = None
								try:
									config = await self.readconfig(job_awaits.path + "config.json")
								except Exception as err:
									print(str(err))
									await job_awaits.Error()

								if config is None:
									await job_awaits.Error()
								elif job_awaits.uid not in self.job_busy:
									await job_awaits.Start()
									await node.jobnode.Transfer(job_awaits)
									if job_awaits.uid not in self.job_busy:
										self.job_busy.add(job_awaits.uid)
										request = {
											'job': job_awaits.__dict__,
											'config': config
										}
										request = Request("execute", **request)
										await websocket.send(request.toJSON())
				message = await self.wsconsumer(node)
				if message is not None:
					request = json.loads(message)
					if "command" in request.keys():
						if request["command"] == "handshake":
							await node.set_node(request["args"]["name"])
							await self.register(node)
							print("Node has been registered: " + node.jobnode.name)
							request = Request("register", **node.jobnode.__dict__)
							await websocket.send(request.toJSON())
						elif request["command"] == "nodestatus":
							node.jobnode.FromData(**request["args"])
							await Database.JobNodeUpdate(node.jobnode)
							if node.jobnode.job_uid > 0 and node.jobnode.progress > 0:
								job = await Database.JobLoad(node.jobnode.job_uid)
								job.extention = node.jobnode.extention
								await job.Progress(node.jobnode.progress)
								print("Job in progress: " + str(node.jobnode.progress))
						elif request["command"] == "jobresult":
							node.jobnode.FromData(**request["args"])
							job = await Database.JobLoad(node.jobnode.job_uid)
							await job.Finish()
							print("Job has been finished: " + str(job.uid))
						else:
							print("unsupported event")
		except (
			websockets.exceptions.ConnectionClosedError,
			websockets.exceptions.ConnectionClosedOK
		) as e:
			await self.unregister(websocket)
			print("Node has been disconnected: " + node.jobnode.name)

			if node.jobnode is not None:
				if node.jobnode.job_uid != 0:
					job = await Database.JobLoad(node.jobnode.job_uid)
					await job.Cancel()

				await node.jobnode.Disconnect()

	async def dbconnect(self):
		try:
			if not os.path.exists(path_root):
				os.mkdir(path_root)
		except Exception as e:
			print('WARNING -- ' + str(e))

		#os.remove(path_root + dbname)
		await Database.Init(os.environ.get('DATABASE_HOST', 'postgres'), os.environ.get('DATABASE_PORT', 5432))
		#await Database.Init('localhost', 5434)
		await Database.Restart()

try:
	astroServer = AstroServer()
	serve_connection = websockets.serve(astroServer.connect, "*", 8080, max_size=1048576*500)
	loop = asyncio.get_event_loop()
	loop.create_task(astroServer.dbconnect())
	loop.run_until_complete(serve_connection)
	print("Server is ready.")
	loop.run_forever()
	#loop.close()
except KeyboardInterrupt:
	print("Server is shuttong down...")
