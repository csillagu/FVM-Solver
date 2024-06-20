import numpy as np
import fluid_fvm.geometry as geo
class MeshConfig():
    pass


class RectangularConfig(MeshConfig):

    # y+1, x-1       y+1  , x             y+1, x+1,
    #       
    #
    #  y, x-1     	y,x      	    	 y, x+1,
    #
    #
    #  y-1,x-1  	y-1,x                y-1, x+1,
    #

    def __init__(self, yNum, xNum, snapLines) -> None:
        self.fxNum =xNum
        self.fyNum = yNum
        self.vxNum =xNum-1
        self.vyNum = yNum-1
        self.snapLines = snapLines
        self.faceMesh = np.zeros((self.fyNum ,self.fxNum),dtype=MeshPoint,)
        self.volumeMesh = np.zeros((self.vyNum, self.vxNum),dtype=MeshPoint,)
        self.connections = []

    def plotMesh(self, ax, vTexts = False, fTexts = False):
        for iy, ix in np.ndindex(self.faceMesh.shape):
            mp = self.faceMesh[iy, ix]
            if isinstance(mp, MeshPoint):
                if fTexts:
                    mp.plot(ax, text="i:"+str(iy)+" j:"+str(ix)+" mid:"+str(self.geo2mathFace((iy,ix))))
                else:
                    mp.plot(ax)
        

        for iy, ix in np.ndindex(self.volumeMesh.shape):
            mp = self.volumeMesh[iy, ix]
            if isinstance(mp, MeshPoint):
                if vTexts:
                    mp.plot(ax, fmt="gx", text="i:"+str(iy)+" j:"+str(ix)+" mid:"+str(self.geo2mathVolume((iy,ix))))
                else:
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
        for k in range(self.fxNum):
            self.faceMesh[0,k] = lineCut(base.lines[0], k/(self.fxNum-1))
        for j in range(self.fyNum):
            for i in range(self.fxNum):
                x = (base.lines[1].p2.x-base.lines[1].p1.x)*j/(self.fyNum-1)
                y = (base.lines[1].p2.y-base.lines[1].p1.y)*j/(self.fyNum-1)
                
                self.faceMesh[j,i] = moveLine(self.faceMesh[0,i], x, y)

    def constructVolumeMesh(self):
        for iy, ix in np.ndindex(self.volumeMesh.shape):
            volPoint = geo.Polygon([geo.Vector(self.faceMesh[iy,ix].x, self.faceMesh[iy,ix].y),
                                    geo.Vector(self.faceMesh[iy+1,ix].x, self.faceMesh[iy+1,ix].y),
                                    geo.Vector(self.faceMesh[iy+1,ix+1].x, self.faceMesh[iy+1,ix+1].y),
                                    geo.Vector(self.faceMesh[iy,ix+1].x, self.faceMesh[iy,ix+1].y)])
            
            self.volumeMesh[iy,ix] = MeshPoint(volPoint.centerMass().x, volPoint.centerMass().y)

    def geo2mathVolume(self, idx):
        # y,x
        y = idx[0]
        x = idx[1]
        return y*(self.vxNum)+x
    
    def math2geoVolume(self, mathidx):
        return (int(np.floor(mathidx/(self.vxNum))), mathidx-int(np.floor(mathidx/(self.vxNum)))*(self.vxNum))
        
    def getVolumeNodeNum(self):
        return self.vyNum*self.vxNum
    
    def getVNode(self, mathid):
        return self.volumeMesh[self.math2geoVolume(mathid)]
    
    def isValidVGeoIdx(self, geoIdx):
        iy = geoIdx[0]
        ix = geoIdx[1]
        return iy>=0 and ix>=0 and iy<self.vyNum and ix<self.vxNum
    
    def getNeighbouringVolumes(self, mathIdx):
        assert(mathIdx>=0)
        assert(mathIdx<=self.getVolumeNodeNum())
        geoIdx = self.math2geoVolume(mathIdx)

        iy = geoIdx[0]
        ix = geoIdx[1]
        ret = []
        for iy_diff, ix_diff in [(-1,0),(0,1),(1,0),(0,-1)]:
                if not self.isValidVGeoIdx((iy+iy_diff,ix+ix_diff)):
                    ret.append([])
                    continue
                ret.append(self.geo2mathVolume((iy+iy_diff,ix+ix_diff)))

        return ret
    
    # Face functions
    def geo2mathFace(self, idx):
        # y,x
        y = idx[0]
        x = idx[1]
        return y*(self.fxNum)+x
    
    def math2geoFace(self, mathidx):
        return (int(np.floor(mathidx/(self.fxNum))), mathidx-int(np.floor(mathidx/(self.fxNum)))*(self.fxNum))
        
    def getFaceNodeNum(self):
        return self.fyNum*self.fxNum
    
    def getFNode(self, mathid):
        return self.faceMesh[self.math2geoFace(mathid)]
    
    def isValidFGeoIdx(self, geoIdx):
        iy = geoIdx[0]
        ix = geoIdx[1]
        return iy>=0 and ix>=0 and iy<self.fyNum and ix<self.fxNum    
    
    def getNeighbouringFaces(self, mathVIdx):
        assert(mathVIdx>=0)
        assert(mathVIdx<=self.getVolumeNodeNum())

        geoIdx = self.math2geoVolume(mathVIdx)

        iy = geoIdx[0]
        ix = geoIdx[1]
        facePointList = []

        for iy_diff, ix_diff in [(0,0), (0,1), (1,1), (1,0)]:
                if not self.isValidFGeoIdx((iy+iy_diff,ix+ix_diff)):
                    raise ValueError("Invalid face mesh found")
                
                facePointList.append(self.geo2mathFace((iy+iy_diff,ix+ix_diff)))

        return facePointList
    
    def getNeighbouringFaceLines(self, mathVIdx):
        facePointList = self.getNeighbouringFaces(mathVIdx=mathVIdx)
        ret = []
        for k in range(len(facePointList)):
            f_line = self.getLineFromFNodes((facePointList[k], facePointList[(k+1)%len(facePointList)]))
            ret.append(f_line)

        return ret
    
    def getLineFromFNodes(self, fNodes):
        mathFId1, mathFId2 = fNodes
        return geo.Line(self.getFNode(mathFId1), self.getFNode(mathFId2))
    
    def getAreaOfElement(self, mathVIdx):
        facePointList = self.getNeighbouringFaces(mathVIdx=mathVIdx)
        x = [self.getFNode(fp).x for fp in facePointList]
        y = [self.getFNode(fp).y for fp in facePointList]
        return 0.5*np.abs(np.dot(x,np.roll(y,1))-np.dot(y,np.roll(x,1)))
        
        



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
    
    def plot(self, ax, fmt = "bx", text = ""):
        ax.plot(self.x, self.y,fmt)
        ax.text( self.x,  self.y, text)

