#!/bin/sh

server="$PWD/modeling_module/physical_problems"
web="$PWD/configurator/configs"
req="$PWD/modeling_module/requirements"

if [ ! -d "$web" ]; then
	mkdir "$web"
fi
rm -rf $web/*

if [ ! -d "$req" ]; then
	mkdir "$req"
fi
rm -rf $req/*


for problem in $(ls "$server"); do
	if [ -f "$server/$problem/config.json" ]; then
		ln -sv "$server/$problem/config.json" "$web/$problem.json"
	fi
	if [ -f "$server/$problem/requirements.txt" ]; then
		cp -v "$server/$problem/requirements.txt" "$req/$problem.txt"
	fi
done
