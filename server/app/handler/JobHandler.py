from database.Database import Database

class JobHandler:

    @staticmethod
    async def Create(session, priority, name):
        user = await Database.UserLoadBySession(session)
        if user is None:
            raise Exception('User not found')
        if not user.IsAuthorized():
            raise Exception('User not authorized')

        path_root = "/var/astromodel/data/"
        job = await Database.JobCreate(user.uid, name)
        path = path_root + str(user.uid) + "/" + str(job.uid) + "/"
        await job.Prepare(path, priority)
        return job

    @staticmethod
    async def GetAll(session):
        user = await Database.UserLoadBySession(session)
        if user is None:
            raise Exception('User not found')
        if not user.IsAuthorized():
            raise Exception('User not authorized')
        return await Database.JobLoadAllByUser(user.uid)

    @staticmethod
    async def Get(session, job_uid):
        user = await Database.UserLoadBySession(session)
        if user is None:
            raise Exception('User not found')
        if not user.IsAuthorized():
            raise Exception('User not authorized')
        return await Database.JobLoad(job_uid)

    @staticmethod
    async def Remove(job):
        return await Database.JobRemove(job.uid)
