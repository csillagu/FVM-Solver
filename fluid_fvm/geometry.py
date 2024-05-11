class Assembly():
    def __init__(self) -> None:
        pass

class Polygon():
    def __init__(self, point_list) -> None:
        self.points = point_list
        
    def getPointsXY(self):
        points_y = [point.y for point in self.points]
        points_x = [point.x for point in self.points]
        return points_x, points_y
        
    def plot(self, ax, fmt="b"):
        points_x, points_y = self.getPointsXY()

        ax.fill(points_x, points_y, color = "#c8c8ff")
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