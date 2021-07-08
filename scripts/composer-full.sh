#!/bin/bash

if [[ ./scripts/composer-full.sh != $BASH_SOURCE ]]; then
  echo "Oh no... You should run this script from repository root: ./scripts/composer-full.sh [ARGS]"
  exit 1
fi


network_name=astro-net
db_filename=astromodel.db
pgdata_host_path=/var/astromodel/pgdata

problems=$PWD/modeling_module/physical_problems
configs=$PWD/configurator/configs
requirements=$PWD/modeling_module/requirements
setups=$PWD/modeling_module/setups


for trashdir in $configs $requirements $setups; do
  echo Creating $trashdir
  if [ ! -d $trashdir ]; then
    mkdir $trashdir
  fi
  rm -rf $trashdir/*
done

for problem in $(ls $problems); do
  if [ -f $problems/$problem/config.json ]; then
    ln -s $problems/$problem/config.json $configs/$problem.json
  fi
  if [ -f $problems/$problem/requirements.txt ]; then
    cp $problems/$problem/requirements.txt $requirements/$problem.txt
  fi
  if [ -f $problems/$problem/setup.sh ]; then
    cp $problems/$problem/setup.sh $setups/$problem
  fi
done


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

# docker system prune -f

docker-compose \
  -p astromodel \
  -f ./docker/postgres/docker-compose.yml \
  -f ./docker/astromodel/docker-compose.yml \
  -f ./docker/wsserver/docker-compose.yml \
  -f ./docker/wsclient/docker-compose.yml \
  -f docker-compose.yml \
  "$@"