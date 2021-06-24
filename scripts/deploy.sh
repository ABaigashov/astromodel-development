#!/usr/bin/env bash

network_name="astro-net"
db_filename="astromodel.db"
pgdata_host_path="/var/astromodel/pgdata"

docker network ls | grep [[:blank:]]$network_name[[:blank:]]

if [[ 1 == $? ]]; then
  echo "Creating network $network_name ..."
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
