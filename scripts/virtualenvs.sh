#!/bin/bash

if [ "$1" == "" ]; then
	for requirement in $(ls ./requirements); do
		venvname=${requirement%%.*}
		python3 -m venv ${PWD}/enviroments/$venvname --system-site-packages
		source ${PWD}/enviroments/$venvname/bin/activate
		pip3 install -r ./requirements/$requirement
		deactivate
	done
else
	pip3 install -r ./requirements/$1.txt
fi