<p align="center">
  <a href="https://astromodel.ru"><img alt="Astromodel" src="./configurator/static/images/logo.svg"></a>
  <a href="https://astromodel.ru"><img alt="Astromodel" src="./configurator/static/images/fond.svg"></a>
</p>

***AstroModel*** project - an astronomical development environment created for automatic storage, deployment, processing and modification of different physical,
astronomical and astrophysical problems containerized by the `Docker`. The project is focused on maximum simplification and
automatization of user interaction with program code, covering a wide range of tasks - numerical modeling,
statistical analysis, optical instrument control, etc.
Any ***AstroModel*** user can easily use,
adapt and reformulate any of the existing tasks or suggest
his own task without the long process of setting up and configuring the software
provision. This allows ***AstroModel*** to act as a versatile
an environment for exchanging proprietary software products in the field of astrophysical computing.


## Installation of Docker and docker-compose on Ubuntu ##

```bash
$ sudo apt update
$ sudo apt upgrade -y
$ sudo apt install curl -y
$ curl -fsSL https://get.docker.com -o get-docker.sh
$ sudo sh get-docker.sh
$ rm get-docker.sh
$ sudo curl -L "https://github.com/docker/compose/releases/download/1.29.2/docker-compose-$(uname -s)-$(uname -m)" -o /usr/local/bin/docker-compose
$ sudo chmod +x /usr/local/bin/docker-compose
$ docker --version && docker-compose --version
```
## Installation of Docker and docker-compose on Mac ##
To install Docker and docker-compose on MacOS, you need<br>
to go [this page](https://docs.docker.com/docker-for-mac/install/)
, choice intel/apple processor version, download and install<br>
docker as normal Mac application. After that, you need to run<br>
this app, close all docker windows and keep this in tray<br>
the icon of docker will be shown in topbar menu<br>

## Some code tips ##

```bash
# If you need to test some code, but you have no generated
# problem parameters you can run "tester.sh". This utility
# is created to run some files with special libraries you have
# in your problem. the testin directory must have "__init__.py" file
# to run tester correctly. In this file you can just import files
# you need like "import file_1" and "file_1.py" is near to init file.

# Run local tester:
$ sudo ./scripts/boot/tester.sh <problem-name> <your-directory>
# For example:
$ sudo ./scripts/boot/tester.sh testing_libs /path/to/leopart/tests
```

```bash
# If you need to create parameters for your problem, just run this
# command. It opens your browser and shows you configuration web page,
# which will be automaticly generated by your "config.yml" file.

# Run web from configuration:
$ sudo ./scripts/boot/config.sh <problem-name>
# For example:
$ sudo ./scripts/boot/config.sh particle_simulator
```

```bash
# If you need to test your config without booting full server
# you can use this utility. Solver can run your problem code
# to test some configs or solve single problem.

# Run local problem solver:
$ sudo ./scripts/boot/solver.sh <problem-name> <parameters.json>
# For example:
$ sudo ./scripts/boot/solver.sh particle_simulator 3d.json
```

```bash
# If you need to have a console of your porblem, like typing
# in command line "python" or "python3" you need to use this
# utility. Command opens interractive console with loaded
# environment and libraries you wrote in your configuration files

# Run interractive console mode:
$ sudo ./scripts/boot/console.sh <problem-name>
# For example:
$ sudo ./scripts/boot/console.sh testing_libs
# To close interractive console mode, write in terminal 
# >>> exit()
```

```bash
# Run server:
$ sudo ./scripts/boot/server.sh [args]
# Acshelly we use this:
$ sudo ./scripts/boot/server.sh up --build

# ATTENTION !!!!
# If you are running server on your computer
# you need to close it CORRECTLY. When you need
# to stop server, press 'CTRL' + 'C' ONLY ONSE.
# (or 'control' + 'C' on MAC). This will take some
# time, but it's ok. If you press this keys one,
# Docker will close faster, but INCORRECTLY.
# !!! PRESS THIS KEYS COMBINATION ONES AND WAIT !!!
# If you accidentally press this keys twice, you
# need to run server again and close it CORRECTLY.
```

```bash
# Clean useless containers:
$ sudo docker system prune -f
$ sudo docker container kill $(docker ps -q)

# Clean all containers:
$ sudo docker system prune -af
```


-----------------------------------------

Compatibility of our project

| **Operating system** | **Compatibility** |
|:--------------------:|:-----------------:|
|    Ubuntu / Debian   |          ✔       |
|      Arch linux      |          ✔       |
|         Mint         |          ✔       |
|        Fedora        |          ✔       |
|         MacOS        |          ✔       |
|      Windows NT      |          ❌      |
|       ChromeOS       |          ❔       |

-----------------------------------------

> Toxicity doesn't lose games, <br>
> &nbsp;&nbsp;&nbsp;bad players do. <br>
> &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;\- Bruno "BCko" Bančić
