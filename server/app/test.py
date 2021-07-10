from database.Database import Database
from handler.JobHandler import JobHandler
from handler.UserHandler import UserHandler
from handler.JobNodeHandler import JobNodeHandler
import os, sys, asyncio

path_root = "/var/astromodel/"
dbname = "astromodel.db"

@asyncio.coroutine
async def test():
    await Database.Init(5432)

    try:
        user = await UserHandler.Authorize('fail_uniq_key', 'developer', 3600)
        print(user)
    except Exception as e:
        print('OK -- ' + str(e))

    try:
        job = await JobHandler.Create('fail_session', 0)
        print(job)
    except Exception as e:
        print('OK -- ' + str(e))

    try:
        user = await UserHandler.Register('dev_uniq_key', 'developer_old', 'Developer Old')
        print('OK -- ' + str(user))
    except Exception as e:
        print('FAIL --' + str(e))

    try:
        await user.Refresh('developer', 'Developer')
        print('OK -- ' + str(user))
    except Exception as e:
        print('FAIL --' + str(e))

    try:
        user = await UserHandler.Authorize('dev_uniq_key', 'developer_session', 3600)
        print('OK -- ' + str(user))
    except Exception as e:
        print('FAIL -- ' + str(e))

    try:
        job = await JobHandler.Create('developer_session', 0)
        print('OK -- ' + str(job))
    except Exception as e:
        print('FAIL -- ' + str(e))

    try:
        job_path = os.path.join(path_root, str(user.uid), str(job.uid))
        await job.Prepare(job_path + "/", 0)
        print('OK -- ' + str(job))
    except Exception as e:
        print('FAIL -- ' + str(e))

    try:
        job = await JobHandler.Get('developer_session', job.uid)
        print('OK -- ' + str(job))
    except Exception as e:
        print('FAIL -- ' + str(e))

    try:
        jobs = await JobHandler.GetAll('developer_session')
        print('OK -- ' + str(jobs))
    except Exception as e:
        print('FAIL -- ' + str(e))

    try:
        node = await JobNodeHandler.Register('test_node_01', 'Test Node (Corei7)', 0)
        print('OK -- ' + str(node))
    except Exception as e:
        print('FAIL -- ' + str(e))

    try:
        node = await JobNodeHandler.Get('test_node_01')
        print('OK -- ' + str(node))
    except Exception as e:
        print('FAIL -- ' + str(e))

    try:
        node = await JobNodeHandler.GetAll()
        print('OK -- ' + str(node))
    except Exception as e:
        print('FAIL -- ' + str(e))

loop = asyncio.get_event_loop()
loop.run_until_complete(test())
#loop.close()
#loop.run_forever()
