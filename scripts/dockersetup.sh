#!/bin/bash

for setup in $(ls ./setups); do
	bash setup
done
