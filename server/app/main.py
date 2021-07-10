from starlette.applications import Starlette
from starlette.responses import JSONResponse
from starlette.responses import Response
from database.Database import Database
from handler.UserHandler import UserHandler
from handler.JobHandler import JobHandler
from handler.JobNodeHandler import JobNodeHandler
import codecs, json, asyncio, requests, base64, shutil
import os, sys

SESSION_DURATION = 3600
app = Starlette(debug=True)
secret = '369XfjmDMsGRRfJp'
#os.remove(path_root + "astromodel.db")
#db = codecs.open(path_root + "astromodel.db", "w+", "utf-8")
#db.close()
#os.chmod(path_root + "astromodel.db", 0o777)

async def dbconnect():
    path_root = "/var/astromodel/"
    if not os.path.exists(path_root):
        os.mkdir(path_root)

    await Database.Init('postgres', 5432)
    #await Database.Init('localhost', 5434)

async def parse_request(request):
    args = await request.body()
    args = args.decode("utf-8")
    return json.loads(args)

async def request_user(email):
    message = 'gudkovla@gmail.com:astromodelgd1810p88'
    message_bytes = message.encode('ascii')
    base64_bytes = base64.b64encode(message_bytes)
    base64_message = base64_bytes.decode('ascii')

    headers = {'Authorization': 'Basic ' + base64_message}
    result = requests.get(
        url = 'https://astromodel.ru/wp-json/wp/v2/users/?search=' + email,
        headers = headers
    )
    return json.loads(result.text)

@app.route("/api/")
async def homepage(request):
    await dbconnect()
    headers = {
        "Content-Type": "text/html"
    }
    return Response('Is nothing here!', 200, headers)

@app.route("/api/user_register", methods=["POST", "GET"])
async def user_register(request):
    await dbconnect()
    headers = {
        "Access-Control-Allow-Origin": "*",
        "Content-Type": "text/json"
    }
    try:
        args = await parse_request(request)
        if 'data' not in args:
            raise Exception("Request is incorrect!")

        data = args['data']
        if ('login' not in data) or ('name' not in data):
            raise Exception("Request is incorrect!")

        result = await request_user(data['login'])
        if len(result) == 0:
            raise Exception("User not found!")

        user = result[0]
        user = await UserHandler.Register(user['id'], data['login'], data['name'])

        response = user.__dict__
        return Response(json.dumps({"answer": response}), 200, headers)
    except Exception as e:
        return Response(json.dumps({"error": str(e)}), 200, headers)

@app.route("/api/user_authorize", methods=["POST"])
async def user_authorize(request):
    await dbconnect()
    headers = {
        "Access-Control-Allow-Origin": "*",
        "Content-Type": "text/json"
    }
    try:
        args = await parse_request(request)
        if 'data' not in args:
            raise Exception("Request is incorrect!")

        data = args['data']
        if ('login' not in data) or ('session' not in data):
            raise Exception("Request is incorrect!")

        result = await request_user(data['login'])
        if len(result) == 0:
            raise Exception("User not found!")

        user = result[0]
        user = await UserHandler.Authorize(user['id'], data['session'], SESSION_DURATION)
        response = user.__dict__
        return Response(json.dumps({"answer": response}), 200, headers)
    except Exception as e:
        return Response(json.dumps({"error": str(e)}), 200, headers)

@app.route("/api/job_execute", methods=["POST"])
async def job_execute(request):
    await dbconnect()
    headers = {
        "Access-Control-Allow-Origin": "*",
        "Content-Type": "text/json"
    }
    try:
        args = await parse_request(request)
        if 'data' not in args:
            raise Exception("Request is incorrect!")

        data = args['data']
        if ('session' not in data) or ('config' not in data) or ('filename' not in data):
            raise Exception("Request is incorrect!")

        job = await JobHandler.Create(data['session'], 0, data['filename'])
        try:
            if not os.path.exists(job.path):
                os.makedirs(job.path)

            file = codecs.open(job.path + "config.json", "w", "utf-8")
            file.write(json.dumps(data['config']))
            file.close()
        except:
            job.Error()

        response = job.__dict__
        await UserHandler.Refresh(data['session'], SESSION_DURATION)
        return Response(json.dumps({"answer": response}), 200, headers)
    except Exception as e:
        return Response(json.dumps({"error": str(e)}), 200, headers)

@app.route("/api/job_check", methods=["POST"])
async def job_check(request):
    await dbconnect()
    headers = {
        "Access-Control-Allow-Origin": "*",
        "Content-Type": "text/json"
    }
    try:
        args = await parse_request(request)
        if 'data' not in args:
            raise Exception("Request is incorrect!")

        data = args['data']
        if ('session' not in data) or ('job_uniq' not in data):
            raise Exception("Request is incorrect!")

        job = await JobHandler.Get(data['session'], data['job_uniq'])
        response = job.__dict__

        await UserHandler.Refresh(data['session'], SESSION_DURATION)
        return Response(json.dumps({"answer": response}), 200, headers)
    except Exception as e:
        return Response(json.dumps({"error": str(e)}), 200, headers)

@app.route("/api/job_list", methods=["POST"])
async def job_list(request):
    await dbconnect()
    headers = {
        "Access-Control-Allow-Origin": "*",
        "Content-Type": "text/json"
    }
    try:
        args = await parse_request(request)
        if 'data' not in args:
            raise Exception("Request is incorrect!")

        data = args['data']
        if ('session' not in data):
            raise Exception("Request is incorrect!")

        jobs = await JobHandler.GetAll(data['session'])

        response = []
        for job in jobs:
            response.append(job.__dict__)

        await UserHandler.Refresh(data['session'], SESSION_DURATION)
        return Response(json.dumps({"answer": response}), 200, headers)
    except Exception as e:
        return Response(json.dumps({"error": str(e)}), 200, headers)

@app.route("/api/job_remove", methods=["POST"])
async def job_list(request):
    await dbconnect()
    headers = {
        "Access-Control-Allow-Origin": "*",
        "Content-Type": "text/json"
    }
    try:
        args = await parse_request(request)
        if 'data' not in args:
            raise Exception("Request is incorrect!")

        data = args['data']
        if ('session' not in data) or ('job_uniq' not in data):
            raise Exception("Request is incorrect!")

        job = await JobHandler.Get(data['session'], data['job_uniq'])
        await JobHandler.Remove(job)

        if os.path.exists(job.path):
            shutil.rmtree(job.path)

        jobs = await JobHandler.GetAll(data['session'])

        response = []
        for job in jobs:
            response.append(job.__dict__)

        await UserHandler.Refresh(data['session'], SESSION_DURATION)
        return Response(json.dumps({"answer": response}), 200, headers)
    except Exception as e:
        return Response(json.dumps({"error": str(e)}), 200, headers)

@app.route("/api/server_status", methods=["GET"])
async def job_list(request):
    await dbconnect()
    headers = {
        "Access-Control-Allow-Origin": "*",
        "Content-Type": "text/json"
    }
    try:
        if ('secret' not in request.query_params):
            raise Exception("Request is incorrect!")

        if request.query_params['secret'] != secret:
            raise Exception("Secret key is incorrect!")

        nodes = await JobNodeHandler.GetAll()

        response = []
        for node in nodes:
            response.append(node.__dict__)

        return Response(json.dumps({"answer": response}), 200, headers)
    except Exception as e:
        return Response(json.dumps({"error": str(e)}), 200, headers)
