#!/bin/bash

if [[ ./scripts/boot/solver.sh != $BASH_SOURCE ]]; then
	echo "Oh no... You should run this script from repository root: ./scripts/boot/solver.sh <problem-name>"
	exit 1
fi
if [ ! -d ./modeling_module/physical_problems/$1 ]; then
	echo "Oh no... No such problem \"$1\", check the spelling"
	exit 1
fi
if [ ! -f ./modeling_module/physical_problems/$1/init_files/$2 ]; then
	echo "Oh no... No such configuration file \"$2\", check the spelling"
	exit 1
fi

bash ./scripts/file-linker.sh $1

export PROBLEM=$1

docker-compose \
	-p astromodel \
	-f ./docker/boot/solver.yml \
	build \
	--build-arg PROBLEM="$1" \
	--build-arg PARAMETERS="$2"

docker-compose \
	-p astromodel \
	-f ./docker/boot/solver.yml \
	up --no-log-prefix

unset PROBLEM
