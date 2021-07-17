#!/bin/bash

if [[ ./scripts/boot/server.sh != $BASH_SOURCE ]]; then
	echo "Oh no... You should run this script from repository root: ./scripts/boot/server.sh [args]"
	exit 1
fi


docker container prune -f --filter "until=12h"
bash ./scripts/file-linker.sh


network_name=astro-net
db_filename=astromodel.db
pgdata_host_path=/var/astromodel/pgdata

docker network ls | grep [[:blank:]]$network_name[[:blank:]]

if [[ 1 == $? ]]; then
	echo Creating network $network_name ...
	docker network create $network_name
	if [[ 0 == $? ]]; then
		echo "Network created"
	else
		echo "Can't create network $network_name"
	fi
else
	echo "Network exists $network_name"
fi


if [[ -e $pgdata_host_path/$db_filename ]]; then
	echo "file exists"
else
	echo "no db file"
fi

docker-compose \
	-p astromodel \
	-f ./docker/postgres/docker-compose.yml \
	-f ./docker/astromodel/docker-compose.yml \
	-f ./docker/wsserver/docker-compose.yml \
	-f ./docker/wsclient/docker-compose.yml \
	-f ./docker/docker-compose.yml \
	"$@"