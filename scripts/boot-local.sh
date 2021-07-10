#!/bin/bash

if [[ ./scripts/boot-local.sh != $BASH_SOURCE ]]; then
  echo "Oh no... You should run this script from repository root: ./scripts/boot-local.sh [ARGS]"
  exit 1
fi


bash ./scripts/file-linker.sh


docker-compose \
  -p astromodel \
  -f ./docker/postgres/docker-compose.yml \
  -f ./docker/astromodel/docker-compose.yml \
  -f ./docker/wsserver/docker-compose.yml \
  -f ./docker/wsclient/docker-compose.yml \
  -f docker-compose.yml \
  "$@"