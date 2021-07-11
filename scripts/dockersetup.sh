#!/bin/bash

for setup in $(ls ./setups); do
	bash ./setups/$setup
done
