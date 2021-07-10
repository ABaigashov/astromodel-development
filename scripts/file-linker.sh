#!/bin/bash

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
    cp $problems/$problem/config.json $configs/$problem.json
  fi
  if [ -f $problems/$problem/requirements.txt ]; then
    cp $problems/$problem/requirements.txt $requirements/$problem.txt
  fi
  if [ -f $problems/$problem/setup.sh ]; then
    cp $problems/$problem/setup.sh $setups/$problem
  fi
done