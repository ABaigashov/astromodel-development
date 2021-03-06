import time

class Job:
	pool = None

	def FromData(self, **kwargs):
		self.uid = kwargs['uid']
		self.user_uid = kwargs['user_uid']
		self.name = kwargs['name']
		self.priority = kwargs['priority']
		self.path = kwargs['path']
		self.extention = kwargs['extention']
		self.problem = kwargs['problem']
		self.state = kwargs['state']
		self.progress = kwargs['progress']
		self.date_created = kwargs['date_created']
		self.date_started = kwargs['date_started']
		self.date_finished = kwargs['date_finished']

	@staticmethod
	async def Init(pool):
		Job.pool = pool
		async with Job.pool.acquire() as conn:
			async with conn.cursor() as cur:
				await cur.execute('''
					CREATE SEQUENCE IF NOT EXISTS seq_jobs_uid
				''')
				await cur.execute('''
					CREATE TABLE IF NOT EXISTS jobs (
						uid BIGINT PRIMARY KEY DEFAULT nextval('seq_jobs_uid'),
						user_uid BIGINT NOT NULL,
						name VARCHAR(32),
						priority INTEGER,
						path VARCHAR(256),
						state VARCHAR(32),
						progress DOUBLE PRECISION,
						date_created BIGINT,
						date_started BIGINT,
						date_finished BIGINT
					)
				''')
				try:
					await cur.execute('''
						ALTER TABLE jobs ADD COLUMN IF NOT EXISTS extention VARCHAR(8)
					''')
					await cur.execute('''
						ALTER TABLE jobs ADD COLUMN IF NOT EXISTS problem VARCHAR(32)
					''')
				except:
					pass

	@staticmethod
	async def Restart():
		async with Job.pool.acquire() as conn:
			async with conn.cursor() as cur:
				await cur.execute('''
					UPDATE jobs SET
						state='await',
						progress=0,
						date_started=0
					WHERE date_finished=0
				''')

	@staticmethod
	async def Create(user_uid, name):
		job = Job()
		job.user_uid = user_uid
		job.priority = 0
		job.name = name
		job.path = ''
		job.extention = ''
		job.problem = ''
		job.state = 'await'
		job.progress = 0
		job.date_created = int(time.time())
		job.date_started = 0
		job.date_finished = 0

		async with Job.pool.acquire() as conn:
			async with conn.cursor() as cur:
				await cur.execute('''
					INSERT INTO jobs (
						user_uid,
						name,
						priority,
						path,
						extention,
						problem,
						state,
						progress,
						date_created,
						date_started,
						date_finished
					)
					VALUES ('{}','{}',{},'{}','{}','{}','{}',{},{},{},{})
					RETURNING uid
				'''.format
					(
					job.user_uid,
					job.name,
					job.priority,
					job.path,
					job.extention,
					job.problem,
					job.state,
					job.progress,
					job.date_created,
					job.date_started,
					job.date_finished
					)
				)

				result = await cur.fetchone()
				job.uid = result[0]
		return job

	@staticmethod
	async def Update(job):
		async with Job.pool.acquire() as conn:
			async with conn.cursor() as cur:
				await cur.execute('''
					UPDATE jobs SET
						user_uid={},
						priority={},
						path='{}',
						extention='{}',
						problem='{}',
						state='{}',
						progress={},
						date_created={},
						date_started={},
						date_finished={}
					WHERE uid={}
				'''.format
					(
					job.user_uid,
					job.priority,
					job.path,
					job.extention,
					job.problem,
					job.state,
					job.progress,
					job.date_created,
					job.date_started,
					job.date_finished,
					job.uid
					)
				)
		return job

	@staticmethod
	async def Load(uniq):
		async with Job.pool.acquire() as conn:
			async with conn.cursor() as cur:
				await cur.execute('''
					SELECT
						uid,
						user_uid,
						name,
						priority,
						path,
						extention,
						problem,
						state,
						progress,
						date_created,
						date_started,
						date_finished
					 FROM jobs WHERE uid={} LIMIT 1'''
					.format(uniq)
				)
				try:
					result = await cur.fetchone()
					job = Job()
					job.uid, job.user_uid, job.name, \
					job.priority, job.path, job.extention, job.problem, \
					job.state, job.progress, job.date_created, \
					job.date_started, job.date_finished = result
					return job
				except Exception as e:
					return None

	@staticmethod
	async def Remove(uid):
		async with Job.pool.acquire() as conn:
			async with conn.cursor() as cur:
				await cur.execute('''
					DELETE FROM jobs WHERE uid={}'''
					.format(uid)
				)
				try:
					result = await cur.fetchone()
					return result[0]
				except Exception as e:
					return None

	@staticmethod
	async def LoadAwaiting(problem):
		async with Job.pool.acquire() as conn:
			async with conn.cursor() as cur:
				await cur.execute('''
					SELECT jb.* FROM
						(SELECT
							uid,
							user_uid,
							name,
							priority,
							path,
							extention,
							problem,
							state,
							progress,
							date_created,
							date_started,
							date_finished
						 FROM jobs
						 WHERE state='await'
						 AND problem='{0}'
						 AND priority=(SELECT MAX(priority) FROM jobs WHERE state='await' AND problem='{0}')
					 ) jb
					 ORDER BY jb.date_created ASC LIMIT 1'''.format(problem))
				try:
					result = await cur.fetchone()
					job = Job()
					job.uid, job.user_uid, job.name, \
					job.priority, job.path, job.extention, job.problem, \
					job.state, job.progress, job.date_created, \
					job.date_started, job.date_finished = result
					return job
				except Exception as e:
					return None

	@staticmethod
	async def LoadAllByUser(user_uid):
		job_list = list()
		async with Job.pool.acquire() as conn:
			async with conn.cursor() as cur:
				await cur.execute('''
					SELECT
						 uid,
						 user_uid,
						 name,
						 priority,
						 path,
						 extention,
						 problem,
						 state,
						 progress,
						 date_created,
						 date_started,
						 date_finished
					FROM jobs
					WHERE user_uid={}
					ORDER BY date_created DESC
				'''.format(user_uid))
				async for row in cur:
					job = Job()
					job.uid, job.user_uid, job.name, \
					job.priority, job.path, job.extention, job.problem, \
					job.state, job.progress, job.date_created, \
					job.date_started, job.date_finished = row
					job_list.append(job)

		return job_list

	async def Prepare(self, path, priority, problem):
		self.path = path
		self.priority = priority
		self.problem = problem
		await Job.Update(self)

	async def Start(self):
		self.state = 'transfer'
		self.progress = 0
		self.date_started = int(time.time())
		await Job.Update(self)

	async def Progress(self, progress):
		self.state = 'execution'
		self.progress = progress
		await Job.Update(self)

	async def Finish(self):
		self.state = 'done'
		self.progress = 1
		self.date_finished = int(time.time())
		await Job.Update(self)

	async def Cancel(self):
		self.state = 'await'
		self.date_started = 0
		self.date_finished = 0
		await Job.Update(self)

	async def Abort(self):
		self.state = 'aborted'
		self.date_finished = int(time.time())
		await Job.Update(self)

	async def Error(self):
		self.state = 'error'
		self.progress = -1
		self.date_finished = int(time.time())
		await Job.Update(self)
