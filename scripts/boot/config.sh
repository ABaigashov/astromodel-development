#!/bin/bash

if [[ ./scripts/boot/config.sh != $BASH_SOURCE ]]; then
	echo "Oh no... You should run this script from repository root: ./scripts/boot/config.sh <problem-name>"
	exit 1
fi
if [ ! -d ./modeling_module/physical_problems/$1 ]; then
	echo "Oh no... No such problem \"$1\", check the spelling"
	exit 1
fi

bash ./scripts/file-linker.sh $1

export PROBLEM=$1

docker-compose \
	-p astromodel \
	-f ./docker/boot/config.yml \
	build

open http://localhost:8008/ || x-www-browser http://localhost:8008/

docker-compose \
	-p astromodel \
	-f ./docker/boot/config.yml \
	up

unset PROBLEM