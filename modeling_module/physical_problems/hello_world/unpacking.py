from constants import constants


class Parser:
    def __init__(self, data):

            self.data = data
            self.objects = []
            self.G = constants["G"]
            self.k = constants["k"]
            self.frames = constants["frames"]
            self.seconds_in_year = constants["seconds_in_year"]
            self.years = constants["years"]
            self.wall = constants["wall"]



    def objects_parse(self):
        for key, obj in self.data.items():
            #try:
                mass = obj["mass"]
                electric_charge = obj["electric_charge"]
                x = obj["coords"][0]
                y = obj["coords"][1]
                vx = obj["vx"]
                vy = obj["vy"]
                particle = Object(mass, electric_charge, x, y, vx, vy)

                self.objects.append(particle)
            #except Exception as e:
               # raise DataError from e
class Object:
    def __init__(self, mass: float, electric_charge: float, x: float, y: float, vx: float, vy: float):
        self.mass = mass
        self.electric_charge = electric_charge
        self.x = x
        self.y = y
        self.vx = vx
        self.vy = vy