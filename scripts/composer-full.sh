#!/bin/bash

if [[ ./scripts/composer-full.sh != $BASH_SOURCE ]]; then
  echo "Oh no... You should run this script from repository root: ./scripts/composer-full.sh [ARGS]"
  exit 1
fi

./scripts/boot.sh

# docker system prune -f

docker-compose \
  -p astromodel \
  -f ./docker/postgres/docker-compose.yml \
  -f ./docker/astromodel/docker-compose.yml \
  -f ./docker/wsserver/docker-compose.yml \
  -f ./docker/wsclient/docker-compose.yml \
  -f docker-compose.yml \
  "$@"