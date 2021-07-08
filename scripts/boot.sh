#!/bin/bash

problems=$PWD/modeling_module/physical_problems
configurations=$PWD/configurator/configs
requirements=$PWD/modeling_module/requirements
network_name=astro-net
db_filename=astromodel.db
pgdata_host_path=/var/astromodel/pgdata

docker network ls | grep [[:blank:]]$network_name[[:blank:]]

if [ ! -d $configurations ]; then
	mkdir $configurations
fi
rm -rf $configurations/*

if [ ! -d $requirements ]; then
	mkdir $requirements
fi
rm -rf $requirements/*


for problem in $(ls $problems); do
	if [ -f $problems/$problem/config.json ]; then
		ln -sv $problems/$problem/config.json $configurations/$problem.json
	fi
	if [ -f $problems/$problem/requirements.txt ]; then
		cp -v $problems/$problem/requirements.txt $requirements/$problem.txt
	fi
	if [ -f $problems/$problemsetup.sh ]; then
		cp -v $problems/$problemsetup.sh $requirements/$problem
	fi
done


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
