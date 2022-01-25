FROM tiangolo/uvicorn-gunicorn-starlette:python3.7
WORKDIR /app

RUN pip install --no-cache-dir aiopg psycopg2 asyncio requests

COPY ${PWD}/server/app /app