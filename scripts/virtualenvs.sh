#!/bin/bash
export WORKON_HOME=$pwd/enviroments
mkdir -p $WORKON_HOME

source /usr/local/bin/virtualenvwrapper.sh

for requirement in $(ls "./requirements"); do
	mkvirtualenv "${requirement%%.*}" -r ./requirements/$requirement
done

mkvirtualenv wsclient -i websockets -i sympy -i numpy -i blessed -i virtualenvwrapper
