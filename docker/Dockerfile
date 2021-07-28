
FROM ubuntu:20.04 as base
WORKDIR /usr/src/

ENV DEBIAN_FRONTEND=noninteractive
RUN apt update && apt upgrade -y

RUN apt-get install python3-dev python3-pip python3-wheel -y
RUN pip3 install --upgrade pip
RUN pip3 install websockets sympy numpy pyyaml

ARG PROBLEM
ENV PROBLEM=$PROBLEM

COPY ${PWD}/modeling_module/setups/$PROBLEM.sh setup.sh
RUN bash setup.sh

COPY ${PWD}/modeling_module/requirements/$PROBLEM.txt requirements.txt
RUN pip3 install -r requirements.txt


FROM base as console

CMD sleep infinity


FROM base as solver

CMD cd folder && python3 -u -B __init__.py


FROM base as problem

COPY ${PWD}/modeling_module/physical_problems/$PROBLEM ./modeling_module/physical_problems/$PROBLEM
COPY ${PWD}/modeling_module/configurator.py ./modeling_module/


FROM problem as client

COPY ${PWD}/modeling_module/jobexecutor.py ./modeling_module/
COPY ${PWD}/server/app ./app

CMD python3 -u -B ./app/wsclient.py
