<p align="center"><a href="https://astromodel.ru"><img alt="Astromodel" src="./configurator/static/images/logo.svg"></a></p>

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
# Run local solver:
$ sudo ./scripts/boot-local.sh <problem-mame> <init file>.json
# For example:
$ sudo ./scripts/boot-local.sh particle_simulator 3d.json


# Run web from configuration:
$ sudo ./scripts/boot-config.sh <problem-mame>
# For example:
$ sudo ./scripts/boot-config.sh particle_simulator


# Run interractive console mode:
$ sudo ./scripts/boot-console.sh <problem-mame>
# For example:
$ sudo ./scripts/boot_console.sh particle_simulator


# Run server:
$ sudo ./scripts/boot-server.sh [args]
# Acshelly we use this:
$ sudo ./scripts/boot-server.sh up --build

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

# Clean useless containers:
$ [sudo] docker system prune -f

# Clean all containers:
$ [sudo] docker system prune -af
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
|      Windows NT      |          ❌       |
|       ChromeOS       |          ❔       |

-----------------------------------------

> Toxicity doesn't lose games, <br>
> &nbsp;&nbsp;&nbsp;bad players do. <br>
> &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;\- Bruno "BCko" Bančić
