FROM python:3

WORKDIR /usr/src/app

RUN pip3 install --no-cache-dir aiopg
RUN pip3 install --no-cache-dir psycopg2
RUN pip3 install --no-cache-dir websockets

CMD [ "python3", "-u", "-B", "wsserver.py" ]
