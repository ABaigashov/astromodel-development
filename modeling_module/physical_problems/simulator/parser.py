from phisycal_constants import contants

class Parser:
    def __init__(self, init):
        self.init = init

    def objects(self):
        for obj in self.init:
            color = self.colors(obj['color'])

    def colors(self):
        pass
