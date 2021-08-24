#!/bin/bash

if [[ ./scripts/boot/tester.sh != $BASH_SOURCE ]]; then
	echo "Oh no... You should run this script from repository root: ./scripts/boot/tester.sh <problem-name> <init-file>.json"
	exit 1
fi
if [ ! -d ./modeling_module/physical_problems/$1 ]; then
	echo "Oh no... No such problem \"$1\", check the spelling"
	exit 1
fi
if [ ! -d $2 ]; then
	echo "Oh no... Directory \"$2\" does not exists, check the spelling"
	exit 1
fi

bash ./scripts/file-linker.sh $1

export FOLDER=$(python -c "import os; print(os.path.realpath('$2'))")

docker-compose \
	-p astromodel \
	-f ./docker/boot/tester.yml \
	build \
	--build-arg PROBLEM="$1"

docker-compose \
	-p astromodel \
	-f ./docker/boot/tester.yml \
	up --no-log-prefix

unset FOLDER
