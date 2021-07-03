#!/bin/bash
for requirement in $(ls ./requirements); do
	venvname=${requirement%%.*}
	python3 -m venv $pwd/enviroments/$venvname
	source $pwd/enviroments/$venvname/bin/activate
	pip3 install -r ./requirements/$requirement
	deactivate
done
