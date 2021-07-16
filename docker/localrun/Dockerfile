FROM ubuntu:20.04
WORKDIR /usr/src/

RUN apt update
RUN apt upgrade -y
ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get install python3-dev python3-pip -y
RUN pip3 install --upgrade pip
RUN pip3 install pyyaml

ARG problem_name
ENV problem_name=$problem_name
ENV problem_dir=./modeling_module/physical_problems/$problem_name

COPY ${PWD}/modeling_module/setups/$problem_name.sh ./setups/$problem_name.sh
COPY ${PWD}/scripts/dockersetup.sh .
CMD ./dockersetup.sh

COPY ${PWD}/modeling_module/requirements/$problem_name.txt ./requirements/$problem_name.txt
COPY ${PWD}/scripts/virtualenvs.sh .
RUN ./virtualenvs.sh $problem_name

COPY ${PWD}/modeling_module ./modeling_module

ARG configuration
ENV configuration=$configuration
CMD python3 -u -B $problem_dir/main.py $(python3 -c "import json,sys,pickle;from modeling_module import configurator;print(pickle.dumps(configurator.Configurator.parse_parameters(json.load(open(sys.argv[1],'rb')))[1]).hex())", $problem_dir/init_files/$configuration) $problem_dir/results result:

