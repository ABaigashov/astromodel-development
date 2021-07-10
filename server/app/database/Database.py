import asyncio, aiopg
from database.Job import Job
from database.User import User
from database.JobNode import JobNode

pg_auth = 'dbname=astromodel user=astromodel password=astromodel host={} port={}'
class Database:
    pool = None

    @staticmethod
    async def Init(host, port):
        if Database.pool is None:
            Database.pool = await aiopg.create_pool(pg_auth.format(host, port))
            await Job.Init(Database.pool)
            await User.Init(Database.pool)
            await JobNode.Init(Database.pool)

    async def Restart():
        await Job.Restart()
        await JobNode.Restart()

    @staticmethod
    async def JobCreate(user_uid, name):
        return await Job.Create(user_uid, name)

    @staticmethod
    async def JobRemove(uid):
        return await Job.Remove(uid)

    @staticmethod
    async def JobLoad(uid):
        return await Job.Load(uid)

    @staticmethod
    async def JobUpdate(job):
        return await Job.Update(job)

    @staticmethod
    async def JobLoadAllByUser(user_uid):
        return await Job.LoadAllByUser(user_uid)

    @staticmethod
    async def JobLoadAwaiting():
        return await Job.LoadAwaiting()

    @staticmethod
    async def UserCreate(external_uid, login, name):
        return await User.Create(external_uid, login, name)

    @staticmethod
    async def UserLoad(uid):
        return await User.Load(uid)

    @staticmethod
    async def UserLoadByExternalUid(external_uid):
        return await User.LoadByExternalUid(external_uid)

    @staticmethod
    async def UserLoadByLogin(login):
        return await User.LoadByLogin(login)

    @staticmethod
    async def UserLoadBySession(session):
        return await User.LoadBySession(session)

    @staticmethod
    async def JobNodeCreate(uniq, name, priority):
        return await JobNode.Create(uniq, name, priority)

    @staticmethod
    async def JobNodeUpdate(node):
        return await JobNode.Update(node)

    @staticmethod
    async def JobNodeLoad(uniq):
        return await JobNode.Load(uniq)

    @staticmethod
    async def JobNodeLoadByName(name):
        return await JobNode.LoadByName(name)

    @staticmethod
    async def JobNodeLoadAll():
        return await JobNode.LoadAll()
