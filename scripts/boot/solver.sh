#!/bin/bash

if [[ ./scripts/boot/solver.sh != $BASH_SOURCE ]]; then
	echo "Oh no... You should run this script from repository root: ./scripts/boot/solver.sh <problem-name> <init-file>.json"
	exit 1
fi
if [ ! -d ./modeling_module/physical_problems/$1 ]; then
	echo "Oh no... No such problem \"$1\", check the spelling"
	exit 1
fi
if [ ! -f ./modeling_module/physical_problems/$1/init_files/$2 ]; then
	echo "Oh no... No such configuration \"$2\", check the spelling"
	exit 1
fi

docker container prune -f --filter "until=12h"
bash ./scripts/file-linker.sh $1

export PROBLEM=$1
export USERVOLUME=$(id -u $USER):$(id -g $USER)

docker-compose \
	-p astromodel \
	-f ./docker/localrun/docker-compose.yml \
	build \
	--build-arg problem_name="$1" \
	--build-arg configuration="$2"

unset PROBLEM
unset USERVOLUME

docker-compose \
	-p astromodel \
	-f ./docker/localrun/docker-compose.yml \
	up --no-log-prefix
