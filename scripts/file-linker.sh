#!/bin/bash

problems=$PWD/modeling_module/physical_problems
configs=$PWD/configurator/configs
requirements=$PWD/modeling_module/requirements
setups=$PWD/modeling_module/setups

link(){
	if [ -f $problems/$1/config.json ]; then
		cp -v $problems/$1/config.json $configs/$1.json
	fi
	if [ -f $problems/$1/requirements.txt ]; then
		cp -v $problems/$1/requirements.txt $requirements/$1.txt
	fi
	if [ -f $problems/$1/setup.sh ]; then
		cp -v $problems/$1/setup.sh $setups/$1.sh
	fi
	if [ ! -d $problems/$1/results ]; then
		mkdir $problems/$1/results
	fi
}

for check in $configs $requirements $setups; do
	echo Creating $trashdir
	if [ ! -d $check ]; then
		mkdir $check
	fi
done

if [ "$1" == "" ]; then
	for problem in $(ls $problems); do
		link $problem
	done
else
	link $1
fi