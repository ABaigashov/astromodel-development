
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


FROM base as tester

CMD python3 -u -B -c "import os; os.chdir('./folder'); import __init__"


FROM base as problem

COPY ${PWD}/modeling_module/physical_problems/$PROBLEM ./modeling_module/physical_problems/$PROBLEM
COPY ${PWD}/modeling_module/configurator.py ./modeling_module/


FROM problem as solver

ARG PARAMETERS
ENV PARAMETERS=$PARAMETERS


CMD python3 -u -B -c "import os; from uuid import uuid4; os.chdir('./modeling_module'); from configurator import Configurator; os.chdir('./physical_problems/$PROBLEM'); from __init__ import Model; os.chdir('/usr/src'); cfg=Configurator('/usr/src/modeling_module/physical_problems/$PROBLEM/init_files/$PARAMETERS'); model=Model(); model.init(cfg, os.path.join('/usr/src/modeling_module/physical_problems/$PROBLEM/results', uuid4().hex[:16]), type('fake_job', (), {'progress': 0})); file=model.run(); os.chmod(file, 0o666);"


FROM problem as client

COPY ${PWD}/modeling_module/jobexecutor.py ./modeling_module/
COPY ${PWD}/server/app ./app

CMD python3 -u -B ./app/wsclient.py
