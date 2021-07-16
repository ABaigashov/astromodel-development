<p align="center"><a href="https://astromodel.ru"><img alt="Astromodel" src="./configurator/static/images/logo.svg"></a></p>

This is very long text and i don't know what to write here. <br>
I think it must to be some information about this project there, <br>
but the repo is private and no one can read this <br>

## Installation on empty Ubuntu 20.04.2.0 LTS ##

```bash
$ sudo apt update
$ sudo apt upgrade -y
$ sudo apt isntall curl -y
$ curl -fsSL https://get.docker.com -o get-docker.sh
$ sudo sh get-docker.sh
$ rm get-docker.sh
$ sudo curl -L "https://github.com/docker/compose/releases/download/1.29.2/docker-compose-$(uname -s)-$(uname -m)" -o /usr/local/bin/docker-compose
$ sudo chmod +x /usr/local/bin/docker-compose
$ docker --version && docker-compose --version
```

## Some code tips ##

```bash
#Run local solver:
$ [sudo] ./scripts/boot-local.sh <problem-mame> <init fle>.json
#For example:
$ [sudo] ./scripts/boot-local.sh particle_simulator 3d.json

#Run server:
$ [sudo] ./scripts/boot-server.sh [args]
#Acshelly we use this:
$ [sudo] ./scripts/boot-server.sh up --build

#Clean useless containers:
$ [sudo] docker system prune -f
#Clean all containers:
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
