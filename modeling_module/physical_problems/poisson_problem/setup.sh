#!/bin/bash

# apt install curl git -y
#
# curl https://cmake.org/files/v3.19/cmake-3.19.3-Linux-x86_64.sh -o cmake.sh
# bash cmake.sh --exclude-subdir --skip-license --prefix=/usr/local

apt-get update
apt-get install --no-install-recommends fenics -y

# git clone https://bitbucket.org/jakob_maljaars/leopart/
# cd leopart/source/cpp
# cmake .
# make
# cd ../..
# python3 setup.py install
# cd ..
