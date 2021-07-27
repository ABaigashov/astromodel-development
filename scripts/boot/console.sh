# #!/bin/bash

# if [[ ./scripts/boot/console.sh != $BASH_SOURCE ]]; then
# 	echo "Oh no... You should run this script from repository root: ./scripts/boot/console.sh <problem-name>"
# 	exit 1
# fi
# if [ ! -d ./modeling_module/physical_problems/$1 ]; then
# 	echo "Oh no... No such problem \"$1\", check the spelling"
# 	exit 1
# fi

# docker container prune -f --filter "until=12h"
# bash ./scripts/file-linker.sh $1

# export PROBLEM=$1

# docker-compose \
# 	-p astromodel \
# 	-f ./docker/console/docker-compose.yml \
# 	build \
# 	--build-arg problem_name="$1" 

# unset PROBLEM

# docker-compose \
# 	-p astromodel \
# 	-f ./docker/console/docker-compose.yml \
# 	up -d --no-log-prefix

# docker exec -it console python3

# docker kill $(docker ps -q --filter ancestor=astromodel_console)