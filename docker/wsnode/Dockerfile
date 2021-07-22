
FROM ubuntu:20.04 as init
WORKDIR /usr/src/

ENV DEBIAN_FRONTEND=noninteractive
RUN apt update && apt upgrade -y

RUN apt-get install python3-dev python3-pip python3-wheel -y
RUN pip3 install --upgrade pip
RUN pip3 install websockets sympy numpy pyyaml

ARG PROBLEM
ENV PROBLEM=$PROBLEM

COPY ${PWD}/modeling_module/setups/${PROBLEM}.sh setup.sh
RUN /bin/bash -c setup.sh

COPY ${PWD}/modeling_module/requirements/${PROBLEM}.txt requirements.txt
RUN pip3 install -r requirements.txt


FROM init as console
CMD python3


FROM init as base

COPY ${PWD}/modeling_module/physical_problems/${PROBLEM} ./modeling_module/physical_problems/${PROBLEM}
COPY ${PWD}/modeling_module/configurator.py ./modeling_module/


FROM base as ws_client

COPY ${PWD}/modeling_module/jobexecutor.py ./modeling_module/
COPY ${PWD}/server/app ./app

CMD python3 -u -B ./app/wsclient.py
