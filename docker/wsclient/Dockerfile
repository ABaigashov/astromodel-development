FROM python:3.7-buster
WORKDIR /usr/src/

RUN pip3 install virtualenvwrapper

COPY ${PWD}/modeling_module/requirements ./requirements
COPY ${PWD}/scripts/virtualenvs.sh .
RUN [ "/bin/bash", "-c", "./virtualenvs.sh"]

COPY ${PWD}/server/app ./app
COPY ${PWD}/modeling_module ./modeling_module

CMD [ "/bin/bash", "-c", "export WORKON_HOME=$pwd/enviroments \
 && source /usr/local/bin/virtualenvwrapper.sh \
 && workon wsclient \
 && python3 -u -B ./app/wsclient.py" ]
