#!/bin/bash

problems=$PWD/modeling_module/physical_problems
configs=$PWD/configurator/configs
requirements=$PWD/modeling_module/requirements
setups=$PWD/modeling_module/setups


for check in $configs $requirements $setups; do
	echo Creating $trashdir
	if [ ! -d $check ]; then
		mkdir $check
		chmod 777 $check
	fi
done

for problem in $(ls $problems); do
	if [ -f $problems/$problem/config.json ]; then
		if [ -f $problems/$problem/requirements.txt ]; then
			cp -v $problems/$problem/requirements.txt $requirements/$problem.txt
		else
			echo "" > $requirements/$problem.txt
		fi
		if [ -f $problems/$problem/setup.sh ]; then
			cp -v $problems/$problem/setup.sh $setups/$problem.sh
			chmod 666 $setups/$problem.sh
		else
			echo "" > $setups/$problem.sh
		fi
		if [ ! -d $problems/$problem/results ]; then
			mkdir $problems/$problem/results
			chmod 777 $problems/$problem/results
		fi

		cp -v $problems/$problem/config.yml $configs/$problem.yml
		chmod 666 $configs/$problem.json $requirements/$problem.txt
	fi
done
