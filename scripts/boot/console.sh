#!/bin/bash

if [[ ./scripts/boot/console.sh != $BASH_SOURCE ]]; then
	echo "Oh no... You should run this script from repository root: ./scripts/boot/console.sh <problem-name>"
	exit 1
fi
if [ ! -d ./modeling_module/physical_problems/$1 ]; then
	echo "Oh no... No such problem \"$1\", check the spelling"
	exit 1
fi

bash ./scripts/file-linker.sh $1

export SEED=$RANDOM$RANDOM$RANDOM

docker-compose \
	-p astromodel \
	-f ./docker/boot/console.yml \
	build \
	--build-arg PROBLEM="$1" 

docker-compose \
	-p astromodel \
	-f ./docker/boot/console.yml \
	up -d

docker exec -it console-${SEED} python3

docker kill $(docker ps -q --filter name=console-${SEED}) >> /dev/null

unset SEED