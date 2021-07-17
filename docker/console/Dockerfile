FROM ubuntu:20.04
WORKDIR /usr/src/

RUN apt update
RUN apt upgrade -y
ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get install python3-dev python3-pip -y
RUN pip3 install --upgrade pip

ARG problem_name
ENV problem_name=$problem_name

COPY ${PWD}/modeling_module/setups/$problem_name.sh ./setups/$problem_name.sh
COPY ${PWD}/scripts/dockersetup.sh .
RUN ./dockersetup.sh

COPY ${PWD}/modeling_module/requirements/$problem_name.txt ./requirements/$problem_name.txt
COPY ${PWD}/scripts/virtualenvs.sh .
RUN ./virtualenvs.sh $problem_name

CMD sleep infinity