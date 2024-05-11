class Assembly():
    def __init__(self) -> None:
        pass

class Polygon():
    def __init__(self, line_list) -> None:
        self.lines = line_list
        self.points_x = [line.p1.x for line in line_list]
        self.points_y = [line.p1.y for line in line_list]
        
    def plot(self, ax, fmt="b"):
        pass

class Line():
    def __init__(self, p1, p2) -> None:
        assert(p1!=p2)
        self.p1 = p1
        self.p2 = p2
    
    def plot(self, ax, fmt="b-"):
        ax.plot([self.p1.x, self.p2.x], [self.p1.y, self.p2.y])


class Point():
    def __init__(self, x, y) -> None:
        self.x = x
        self.y = y

    def getPosition(self):
        return (self.x, self.y)
    
    def plot(self, ax, fmt = "bx"):
        ax.plot(self.x, self.y,fmt)