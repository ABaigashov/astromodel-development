<p align="center"><a href="https://astromodel.ru"><img alt="Astromodel" src="./configurator/static/images/logo.svg"></a></p>

# code stuff #

```bash

#Install Docker and docker-compose on Ubuntu:
$ sudo apt update
$ sudo apt isntall curl
$ curl -fsSL https://get.docker.com -o get-docker.sh
$ sudo sh get-docker.sh
$ rm get-docker.sh
$ sudo curl -L "https://github.com/docker/compose/releases/download/1.29.2/docker-compose-$(uname -s)-$(uname -m)" -o /usr/local/bin/docker-compose
$ sudo chmod +x /usr/local/bin/docker-compose
$ docker --version && docker-compose --version


#Run local solver:
$ ./scripts/boot-local.sh <problem-mame> <init fle>.json
#For example:
$ ./scripts/boot-local.sh particle_simulator 3d.json

#Run server:
$ ./scripts/boot-server.sh [args]
#Acshelly we use this:
$ ./scripts/boot-server.sh up --build

#Clean useless containers:
$ docker system prune -f
#Clean all containers:
$ docker system prune -af
```


-----------------------------------------

| **Operating system** | **Compatibility** |
|:--------------------:|:-----------------:|
|    Ubuntu / Debian   |          ✔        |
|      Arch linux      |          ✔        |
|         Mint         |          ✔        |
|        Fedora        |          ✔        |
|         MacOS        |          ✔        |
|      Windows NT      |          ❌        |
|       ChromeOS       |          ❔        |

-----------------------------------------

> Toxicity doesn't lose games, <br>
> bad players do. <br>
> &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;\- Bruno "BCko" Bančić
