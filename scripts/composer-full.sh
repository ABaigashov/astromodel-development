#!/bin/bash

result=${PWD##*/}
selfname=${0##*/}

if [[ "astro-model" != "$result" ]]; then
  echo "Oh no... You should run this script from repository root: ./scripts/$selfname [ARGS]"
  exit 1
fi

#ln -sf ${PWD}/modeling_module/particle_simulator/requirements.txt ./docker/wsclient/requirements.txt
#cp ./modeling_module/particle_simulator/requirements.txt ./docker/wsclient/requirements.txt

./scripts/boot.sh

docker-compose \
  -p astromodel \
  -f ./docker/postgres/docker-compose.yml \
  -f ./docker/astromodel/docker-compose.yml \
  -f ./docker/wsserver/docker-compose.yml \
  -f ./docker/wsclient/docker-compose.yml \
  -f docker-compose.yml \
  "$@"
