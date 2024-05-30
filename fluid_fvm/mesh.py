import numpy as np
import fluid_fvm.geometry as geo
class MeshConfig():
    pass


class RectangularConfig(MeshConfig):
    def __init__(self, xNum, yNum, snapLines) -> None:
        self.iNum =xNum
        self.jNum = yNum
        self.snapLines = snapLines
        self.faceMesh = np.zeros((yNum,xNum),dtype=MeshPoint,)
        self.volumeMesh = np.zeros((yNum-1,xNum-1),dtype=MeshPoint,)
        self.connections = []

    def plotMesh(self, ax):
        for iy, ix in np.ndindex(self.faceMesh.shape):
            mp = self.faceMesh[iy, ix]
            if isinstance(mp, MeshPoint):
                mp.plot(ax)
        

        for iy, ix in np.ndindex(self.volumeMesh.shape):
            mp = self.volumeMesh[iy, ix]
            if isinstance(mp, MeshPoint):
                mp.plot(ax, fmt="gx")
    
    def constructMesh(self, base):
        if len(base.lines) != 4:
            raise ValueError("Base must be rectangular")
        #print(base.lines[0].isPerpendicular(base.lines[1]) and base.lines[1].isPerpendicular(base.lines[2]) and base.lines[2].isPerpendicular(base.lines[3]))
        if not (base.lines[0].isPerpendicular(base.lines[1]) and base.lines[1].isPerpendicular(base.lines[2]) and base.lines[2].isPerpendicular(base.lines[3])):
            raise ValueError("Base must be rectangular")
        self.constructFaceMesh(base)
        self.constructVolumeMesh()

        
    def constructFaceMesh(self, base):
        for k in range(self.iNum):
            self.faceMesh[0,k] = lineCut(base.lines[0], k/(self.iNum-1))
        for j in range(self.jNum):
            for i in range(self.iNum):
                x = (base.lines[1].p2.x-base.lines[1].p1.x)*j/(self.jNum-1)
                y = (base.lines[1].p2.y-base.lines[1].p1.y)*j/(self.jNum-1)
                
                self.faceMesh[j,i] = moveLine(self.faceMesh[0,i], x, y)

    def constructVolumeMesh(self):
        for iy, ix in np.ndindex(self.volumeMesh.shape):
            volPoint = geo.Polygon([geo.Point(self.faceMesh[iy,ix].x, self.faceMesh[iy,ix].y),
                                    geo.Point(self.faceMesh[iy+1,ix].x, self.faceMesh[iy+1,ix].y),
                                    geo.Point(self.faceMesh[iy+1,ix+1].x, self.faceMesh[iy+1,ix+1].y),
                                    geo.Point(self.faceMesh[iy,ix+1].x, self.faceMesh[iy,ix+1].y)])
            
            self.volumeMesh[iy,ix] = MeshPoint(volPoint.centerMass().x, volPoint.centerMass().y)

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