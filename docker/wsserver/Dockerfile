FROM python:3
WORKDIR /app

RUN pip3 install --no-cache-dir aiopg psycopg2 websockets
COPY ${PWD}/server/app /app

CMD [ "python3", "-u", "-B", "wsserver.py" ]
