class Assembly():
    def __init__(self) -> None:
        pass

class Polygon():
    def __init__(self, line_list) -> None:
        self.lines = line_list

class Line():
    def __init__(self, p1, p2) -> None:
        assert(p1!=p2)
        self.p1 = p1
        self.p2 = p2


class Point():
    def __init__(self, x, y) -> None:
        self.x = x
        self.y = y

    def getPosition(self):
        return (self.x, self.y)