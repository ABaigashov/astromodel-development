FROM ubuntu:20.04
WORKDIR /usr/src/

RUN apt update
RUN apt upgrade -y

RUN  apt-get install python3-dev python3-pip python3-venv python3-wheel -y
RUN pip3 install virtualenv websockets sympy numpy


COPY ${PWD}/modeling_module/requirements ./requirements
COPY ${PWD}/scripts/virtualenvs.sh .
RUN ./virtualenvs.sh

#COPY ${PWD}/server/app ./app
#COPY ${PWD}/modeling_module ./modeling_module

#CMD [ "/bin/bash", "-c", "export WORKON_HOME=$pwd/enviroments \
# && source /usr/local/bin/virtualenvwrapper.sh \
# && workon wsclient \
# && python3 -u -B ./app/wsclient.py" ]
