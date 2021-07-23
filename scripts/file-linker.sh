#!/bin/bash

problems=$PWD/modeling_module/physical_problems
configs=$PWD/configurator/configs
requirements=$PWD/modeling_module/requirements
setups=$PWD/modeling_module/setups


for check in $configs $requirements $setups; do
	echo Creating $check
	if [ ! -d $check ]; then
		mkdir $check
		chmod 777 $check
	fi
done

cp -v $PWD/docker/docker-compose.template $PWD/docker/docker-compose.yml
chmod 666 $PWD/docker/docker-compose.yml


for problem in $(ls $problems); do
	if [ -f $problems/$problem/config.yml ]; then
		if [ -f $problems/$problem/requirements.txt ]; then
			cp -v $problems/$problem/requirements.txt $requirements/$problem.txt
		else
			echo "" > $requirements/$problem.txt
		fi
		if [ -f $problems/$problem/setup.sh ]; then
			cp -v $problems/$problem/setup.sh $setups/$problem.sh
			chmod 666 $setups/$problem.sh
		else
			echo "" > $setups/$problem.sh
		fi
		if [ ! -d $problems/$problem/results ]; then
			mkdir $problems/$problem/results
			chmod 777 $problems/$problem/results
		fi

		cp -v $problems/$problem/config.yml $configs/$problem.yml
		chmod 666 $configs/$problem.yml $requirements/$problem.txt

		echo ""                                            >> $PWD/docker/docker-compose.yml
		echo "  wsclient-$problem:"                       >> $PWD/docker/docker-compose.yml
		echo "    depends_on:"                             >> $PWD/docker/docker-compose.yml
		echo "      - wsserver"                            >> $PWD/docker/docker-compose.yml
		echo "    volumes:"                                >> $PWD/docker/docker-compose.yml
		echo "      - /var/astromodel:/var/astromodel"     >> $PWD/docker/docker-compose.yml
		echo "    deploy:"                                 >> $PWD/docker/docker-compose.yml
		echo "      replicas: 4"                           >> $PWD/docker/docker-compose.yml
		echo "      restart_policy:"                       >> $PWD/docker/docker-compose.yml
		echo "        condition: on-failure"               >> $PWD/docker/docker-compose.yml
		echo "    build:"                                  >> $PWD/docker/docker-compose.yml
		echo "      context: \${PWD}"                      >> $PWD/docker/docker-compose.yml
		echo "      dockerfile: \${PWD}/docker/Dockerfile" >> $PWD/docker/docker-compose.yml
		echo "      target: client"                        >> $PWD/docker/docker-compose.yml
		echo "      args:"                                 >> $PWD/docker/docker-compose.yml
		echo "        - PROBLEM=$problem"                  >> $PWD/docker/docker-compose.yml
	fi
done
