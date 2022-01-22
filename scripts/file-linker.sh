#!/bin/bash

problems=$PWD/modeling_module/physical_problems
configs=$PWD/frontend/configs
requirements=$PWD/modeling_module/requirements
setups=$PWD/modeling_module/setups
server=$PWD/docker/server_config.yml


for check in $configs $requirements $setups; do
	echo Creating $check
	if [ ! -d $check ]; then
		mkdir $check
		chmod 777 $check
	fi
done

cp $PWD/docker/server.template $PWD/docker/docker-compose.yml
chmod 666 $PWD/docker/docker-compose.yml

docker container prune -f --filter "until=12h"

function link {

	if [ -f $problems/$1/config.yml ]; then
		if [ -f $problems/$1/requirements.txt ]; then
			cp $problems/$1/requirements.txt $requirements/$1.txt
		else
			echo "" > $requirements/$1.txt
		fi
		if [ -f $problems/$1/setup.sh ]; then
			cp $problems/$1/setup.sh $setups/$1.sh
		else
			echo "" > $setups/$1.sh
		fi
		if [ ! -d $problems/$1/results ]; then
			mkdir $problems/$1/results
			chmod 777 $problems/$1/results
		fi
	fi

	cp $problems/$1/config.yml $configs/$1.yml
	chmod 666 $configs/$1.yml $requirements/$1.txt $setups/$1.sh
}

if [[ $1 ]]; then

	link $1

else

	tmp_IFS=$IFS
	IFS=$'\n'

	for j in $(cat $server); do

		IFS=": "
		q=($j)
		problem=${q[0]}
		size=${q[1]}

		link $problem

		if [ -f $problems/$problem/config.yml ]; then
			echo ""                                                  >> $PWD/docker/docker-compose.yml
			echo "    wsclient-$problem:"                            >> $PWD/docker/docker-compose.yml
			echo "        depends_on:"                               >> $PWD/docker/docker-compose.yml
			echo "            - wsserver"                            >> $PWD/docker/docker-compose.yml
			echo "        volumes:"                                  >> $PWD/docker/docker-compose.yml
			echo "            - /var/astromodel:/var/astromodel"     >> $PWD/docker/docker-compose.yml
			echo "        deploy:"                                   >> $PWD/docker/docker-compose.yml
			echo "            replicas: $size"                       >> $PWD/docker/docker-compose.yml
			echo "            restart_policy:"                       >> $PWD/docker/docker-compose.yml
			echo "                condition: on-failure"             >> $PWD/docker/docker-compose.yml
			echo "        build:"                                    >> $PWD/docker/docker-compose.yml
			echo "            context: \${PWD}"                      >> $PWD/docker/docker-compose.yml
			echo "            dockerfile: \${PWD}/docker/Dockerfile" >> $PWD/docker/docker-compose.yml
			echo "            target: client"                        >> $PWD/docker/docker-compose.yml
			echo "            args:"                                 >> $PWD/docker/docker-compose.yml
			echo "                - PROBLEM=$problem"                >> $PWD/docker/docker-compose.yml
		fi
	done

	IFS=$tmp_IFS
	unset tmp_IFS
fi
