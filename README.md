<p align="center">
  <a href="https://astromodel.ru">
    <img alt="Astromodel logo" src="./configurator/static/images/logo.svg">
    <img alt="Presidential grand" src="./configurator/static/images/fond.svg">
  </a>
</p>
<br>

------------------------------------------

<br>
**AstroModel** project is astronomical development environment created for automatic storage,
deployment, processing and modification of different physical, astronomical and astrophysical
problems containerized by the `Docker`. The main goal of this project is simplification and
automatization of user interaction with program code. It covers a wide range of calculus such as
numerical modeling, statistical analysis, optical instrument control, etc.

**AstroModel** user can easily
use, adapt and reformulate any of the existing problems. Also user can suggest his own problem solver
without the long process of setting up and configuring the software provision. This allows our project
to be universal environment for exchanging software products in astrophysical computing scope.
<br>

------------------------------------------

<br>


# Compatibility of our project #

| **Operating system** | **Compatibility** |
|:--------------------:|:-----------------:|
|    Ubuntu / Debian   |          ✔        |
|      Arch linux      |          ✔        |
|         Mint         |          ✔        |
|        Fedora        |          ✔        |
|         MacOS        |          ✔        |
|      Windows NT      |         ❌        |
|       ChromeOS       |         ❔        |


<br><br>


# Our code utilities #

<br>
### *test local files using problem enviroment* ###
If you need to test some code, but you have no generated<br>
problem parameters you can run `tester.sh`. This utility<br>
is created to run some files with special libraries you have<br>
in your problem. the testin directory must have `__init__.py` file<br>
to run tester correctly. In this file you can just import files<br>
you need like `import file_1` and `file_1.py` is near to init file.<br>
```bash
# Run local tester:
$ sudo ./scripts/boot/tester.sh <problem-name> <your-directory>
# For example:
$ sudo ./scripts/boot/tester.sh testing_libs /path/to/leopart/tests
```

<br>
### *create configuration file by web page* ###
If you need to create parameters for your problem, just run this<br>
command. It opens your browser and shows you configuration web page,<br>
which will be automaticly generated by your `config.yml` file.<br>
```bash
# Run web from configuration:
$ sudo ./scripts/boot/config.sh <problem-name>
# For example:
$ sudo ./scripts/boot/config.sh particle_simulator
```

<br>
### *test configuration file without server* ###
If you need to test your config without booting full server<br>
you can use this utility. Solver can run your problem code<br>
to test some configs or solve single problem.<br>
```bash
# Run local problem solver:
$ sudo ./scripts/boot/solver.sh <problem-name> <parameters.json>
# For example:
$ sudo ./scripts/boot/solver.sh particle_simulator 3d.json
```

<br>
### *fast check some enviroment tips by console* ###
If you need to have a console of your porblem, like typing<br>
in command line `python` or `python3` you need to use this<br>
utility. Command opens interractive console with loaded<br>
environment and libraries you wrote in your configuration files.<br>
To close interractive console mode, write in terminal `exit()`.<br>
```bash
# Run interractive console mode:
$ sudo ./scripts/boot/console.sh <problem-name>
# For example:
$ sudo ./scripts/boot/console.sh testing_libs
```

<br>
### *run your own local server* ###
> *!!! ATTENTION !!!* <br>
> If you are running server on your computer you need to close it CORRECTLY. When you need<br>
> to stop server, press 'CTRL' + 'C' ONLY ONCE. (or 'control' + 'C' on Mac). This will take<br>
> some time, but it's ok. If you press this keys one more time, Docker will close faster,<br>
> but INCORRECTLY. You need to press this keys ONLY ONCE and WAIT a little bit. If you<br>
> accidentally press this keys twice, you need to run server again and close it CORRECTLY.<br>
If you need to run server like in our site, you need to use<br>
this command. It makes all things he need to be exactly like<br>
on our server.<br>
```bash
# Run server:
$ sudo ./scripts/boot/server.sh [args]
# Acshelly we use this:
$ sudo ./scripts/boot/server.sh up --build
```

<br>
### *other code tips when something goes wrong* ###
```bash
# Look inside the container:
$ sudo docker exec -it <container-name> bash

# Clean useless containers:
$ sudo docker system prune -f
$ sudo docker container kill $(docker ps -q)

# Clean all containers:
$ sudo docker system prune -af
```


<br><br>


# Installation of Docker and docker-compose on Ubuntu #

```bash
$ sudo apt update
$ sudo apt upgrade -y
$ sudo apt install curl -y
$ curl -fsSL https://get.docker.com -o get-docker.sh
$ sudo sh get-docker.sh
$ rm get-docker.sh
$ version="1.29.2/docker-compose-$(uname -s)-$(uname -m)"
$ url="https://github.com/docker/compose/releases/download/$version"
$ sudo curl -L $url -o /usr/local/bin/docker-compose
$ sudo chmod +x /usr/local/bin/docker-compose
$ docker --version && docker-compose --version
```

<br>

# Installation of Docker and docker-compose on Mac #
To install Docker and docker-compose on MacOS, you need to go [this page](https://docs.docker.com/docker-for-mac/install/),
choice intel / apple<br> processor version, download and install docker as typical Mac application. After that, you need
to run<br> this app, close all docker windows and keep this in tray the icon of `Docker` will be shown in topbar menu<br>


<br><br>


-----------------------------------------

> Toxicity doesn't lose games,<br>
> &nbsp;&nbsp;&nbsp;bad players do.<br>
> &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;\- Bruno "BCko" Bančić
