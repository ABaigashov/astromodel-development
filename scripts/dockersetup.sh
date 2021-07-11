#!/bin/bash

if [ "$1" == "" ]; then
	for setup in $(ls ./setups); do
		bash ./setups/$setup
	done
else
	bash ./setups/$1.sh
fi