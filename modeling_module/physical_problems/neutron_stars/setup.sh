#!/bin/bash

# This file has special name 'setup.sh' and
# you SHOULDN'T RENAME it. This file setup file
# and it will be called on the start of problem
# calculation. It's written in bash language
# and contains all installation information
# of some non-pip libraries like 'fenics',
# 'dolfin', 'mshr', 'leopart' and other...

# [Hint]:
#    Look at the first line. There you can see some
#    strange comment. This is special flag to linux
#    system which says that this file is bash script
#    You ALWAYS NEED to have this commnet in FIRST LINE

# [Hint]:
#    This file will be fed to the 'Docker container'
#    which runs as root user. You SHOULDN'T WRITE 'sudo'
#    in this file. Otherwise error  will be raised

echo '"Hello World!" is installing now!'


# nothing more needed
# check "testing_libs" problem