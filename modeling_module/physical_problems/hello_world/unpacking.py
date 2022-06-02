from constants import constants


class Parser:
    def __init__(self, data, logger):
        self.logger = logger
        self.logger.debug('Start unpacking')
        self.data = data
        self.objects = []
        self.num = 0
        try:
            self.G = constants["G"]
            self.k = constants["k"]
            self.frames = constants["frames"]
            self.seconds_in_year = constants["seconds_in_year"]
            self.years = constants["years"]
            self.wall = constants["wall"]
            self.logger.info('Constants are OK')

        except Exception as e:
            self.logger.error(f"Constants have wrong format: {e} ")

    def objects_parse(self):
        try:
            for key, obj in self.data.items():
                self.num += 1
                mass = obj["mass"]
                electric_charge = obj["electric_charge"]
                x = obj["coords"][0]
                y = obj["coords"][1]
                vx = obj["vx"]
                vy = obj["vy"]
                particle = Object(mass, electric_charge, x, y, vx, vy)

                self.objects.append(particle)
                self.logger.info(f"Create object {self.num}")

        except Exception as e:
            self.logger.error(f"Data error, check input parametrs: {e} ")


class Object:
    def __init__(self, mass: float, electric_charge: float, x: float, y: float, vx: float, vy: float):
        self.mass = mass
        self.electric_charge = electric_charge
        self.x = x
        self.y = y
        self.vx = vx
        self.vy = vy
