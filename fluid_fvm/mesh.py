import numpy as np

class MeshConfig():
    pass


class RectangularConfig(MeshConfig):
    def __init__(self, xNum, yNum, snapLines) -> None:
        self.iNum =xNum
        self.jNum = yNum
        self.snapLines = snapLines
        self.meshPoints = np.zeros((yNum,xNum),dtype=MeshPoint,)

        self.connections = []

    def plotMesh(self, ax):
        for iy, ix in np.ndindex(self.meshPoints.shape):
            mp = self.meshPoints[iy, ix]
            if isinstance(mp, MeshPoint):
                mp.plot(ax)
    
    def constructMesh(self, base):
        if len(base.lines) != 4:
            raise ValueError("Base must be rectangular")
        #print(base.lines[0].isPerpendicular(base.lines[1]) and base.lines[1].isPerpendicular(base.lines[2]) and base.lines[2].isPerpendicular(base.lines[3]))
        if not (base.lines[0].isPerpendicular(base.lines[1]) and base.lines[1].isPerpendicular(base.lines[2]) and base.lines[2].isPerpendicular(base.lines[3])):
            raise ValueError("Base must be rectangular")
        

        for k in range(self.iNum):
            self.meshPoints[0,k] = lineCut(base.lines[0], k/(self.iNum-1))
        for j in range(self.jNum):
            for i in range(self.iNum):
                x = (base.lines[1].p2.x-base.lines[1].p1.x)*j/(self.jNum-1)
                y = (base.lines[1].p2.y-base.lines[1].p1.y)*j/(self.jNum-1)
                
                self.meshPoints[j,i] = moveLine(self.meshPoints[0,i], x, y)
            
def moveLine(point, x,y):
    return MeshPoint(point.x+x, point.y+y)

def lineCut(line, proportion):
    diff_vec = line.p2-line.p1
    return MeshPoint(line.p1.x+diff_vec.x*proportion, line.p1.y+diff_vec.y*proportion)

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