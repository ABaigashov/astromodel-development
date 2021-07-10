import time

class JobNode:
    pool = None

    def FromData(self, **kwargs):
        self.uniq = kwargs['uniq']
        self.name = kwargs['name']
        self.state = kwargs['state']
        self.extention = kwargs['extention']
        self.priority = kwargs['priority']
        self.job_uid = kwargs['job_uid']
        self.progress = kwargs['progress']
        self.date_updated = kwargs['date_updated']

    @staticmethod
    async def Init(pool):
        JobNode.pool = pool
        async with JobNode.pool.acquire() as conn:
            async with conn.cursor() as cur:
                await cur.execute('''
                    CREATE SEQUENCE IF NOT EXISTS seq_job_nodes_uid
                ''')
                await cur.execute('''
                    CREATE TABLE IF NOT EXISTS job_nodes (
                        uid BIGINT PRIMARY KEY DEFAULT nextval('seq_job_nodes_uid'),
                        uniq VARCHAR(36) NOT NULL,
                        name VARCHAR(32) NOT NULL,
                        state VARCHAR(32),
                        priority INTEGER,
                        job_uid BIGINT,
                        progress DOUBLE PRECISION,
                        date_updated BIGINT
                    )
                ''')

    @staticmethod
    async def Restart():
        async with JobNode.pool.acquire() as conn:
            async with conn.cursor() as cur:
                await cur.execute("DROP TABLE IF EXISTS job_nodes")
        await JobNode.Init(JobNode.pool)

    @staticmethod
    async def Create(uniq, name, priority):
        timestamp = int(time.time())
        node = JobNode()
        node.uniq = uniq
        node.name = name
        node.state = 'offline'
        node.priority = priority
        node.job_uid = 0
        node.progress = 0
        node.date_updated = timestamp

        async with JobNode.pool.acquire() as conn:
            async with conn.cursor() as cur:
                await cur.execute('''
                    INSERT INTO job_nodes (
                        uniq,
                        name,
                        state,
                        priority,
                        job_uid,
                        progress,
                        date_updated
                    )
                    VALUES ('{}','{}','{}',{},{},{},{})
                    RETURNING uid
                    '''.format
                    (
                        node.uniq,
                        node.name,
                        node.state,
                        node.priority,
                        node.job_uid,
                        node.progress,
                        node.date_updated
                    )
                )
                result = await cur.fetchone()
                node.uid = result[0]
        return node

    @staticmethod
    async def Update(node):
        async with JobNode.pool.acquire() as conn:
            async with conn.cursor() as cur:
                await cur.execute('''
                    UPDATE job_nodes SET
                        name='{}',
                        state='{}',
                        priority={},
                        job_uid={},
                        progress={},
                        date_updated={}
                    WHERE uid={}
                    '''.format
                    (
                    node.name,
                    node.state,
                    node.priority,
                    node.job_uid,
                    node.progress,
                    int(time.time()),
                    node.uid,
                    )
                )
        return node

    @staticmethod
    async def Load(uniq):
        async with JobNode.pool.acquire() as conn:
            async with conn.cursor() as cur:
                await cur.execute('''
                    SELECT * FROM job_nodes WHERE uniq='{}' LIMIT 1
                '''.format(uniq))
                async for row in cur:
                    try:
                        node = JobNode()
                        node.uid = row[0]
                        node.uniq = row[1]
                        node.name = row[2]
                        node.state = row[3]
                        node.priority = row[4]
                        node.job_uid = row[5]
                        node.progress = row[6]
                        node.date_updated = row[7]
                        return node
                    except Exception as e:
                        return None

    @staticmethod
    async def LoadByName(name):
        async with JobNode.pool.acquire() as conn:
            async with conn.cursor() as cur:
                await cur.execute('''
                    SELECT * FROM job_nodes WHERE name='{}' LIMIT 1
                '''.format(name));
                async for row in cur:
                    try:
                        node = JobNode()
                        node.uid = row[0]
                        node.uniq = row[1]
                        node.name = row[2]
                        node.state = row[3]
                        node.priority = row[4]
                        node.job_uid = row[5]
                        node.progress = row[6]
                        node.date_updated = row[7]
                        return node
                    except Exception as e:
                        return None

    @staticmethod
    async def LoadAll():
        node_list = list()
        async with JobNode.pool.acquire() as conn:
            async with conn.cursor() as cur:
                await cur.execute('''
                    SELECT * FROM job_nodes
                    ORDER BY priority,name DESC
                ''')
                async for row in cur:
                    node = JobNode()
                    node.uniq = row[0]
                    node.uniq = row[1]
                    node.name = row[2]
                    node.state = row[3]
                    node.priority = row[4]
                    node.job_uid = row[5]
                    node.progress = row[6]
                    node.date_updated = row[7]
                    node_list.append(node)

        return node_list

    async def Refresh(self, name, priority):
        self.name = name
        self.priority = priority
        self.date_updated = int(time.time())
        await JobNode.Update(self)

    async def Transfer(self, job):
        self.state = 'transfer'
        self.job_uid = job.uid
        self.progress = job.progress
        self.date_updated = int(time.time())
        await JobNode.Update(self)

    async def Disconnect(self):
        self.state = 'offline'
        self.job_uid = 0
        self.date_updated = int(time.time())
        await JobNode.Update(self)
