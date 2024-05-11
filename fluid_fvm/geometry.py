class Assembly():
    def __init__(self, base_polygon, polygon_list) -> None:
        self.polygon_list = polygon_list
        self.base_polygon = base_polygon

    def plot(self, ax):
        self.base_polygon.plot(ax, "white")
        for p in self.polygon_list:
            p.plot(ax)
        

class Polygon():
    def __init__(self, point_list) -> None:
        self.points = point_list
        
    def getPointsXY(self):
        points_y = [point.y for point in self.points]
        points_x = [point.x for point in self.points]
        return points_x, points_y
        
    def plot(self, ax, color="#c8ffc8"):
        points_x, points_y = self.getPointsXY()

        ax.fill(points_x, points_y, color = color, edgecolor="black")
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