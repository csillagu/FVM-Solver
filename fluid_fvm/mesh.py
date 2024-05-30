class MeshConfig():
    pass


class RectangularConfig(MeshConfig):
    def __init__(self, xNum, yNum, snapLines) -> None:
        self.xNum =xNum
        self.yNum = yNum
        self.snapLines = snapLines
        self.meshPoints = []
        self.connections = []

    def plotMesh(self, ax):
        for p in self.meshPoints:
            p.plot(ax)




class MeshPoint():
    def __init__(self, x, y) -> None:
        self.x = x
        self.y = y

    def __eq__(self, other: object) -> bool:
        return (self.x == other.x) and (self.y == other.y)
    
    def __sub__(self, other):
        #self-other
        return self.x-other.x, self.y-other.y
    
    def plot(self, ax, fmt = "bx"):
        ax.plot(self.x, self.y,fmt)