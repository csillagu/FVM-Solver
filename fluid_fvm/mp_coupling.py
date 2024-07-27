import fluid_fvm.physics as ph
import fluid_fvm.geometry as geo

import numpy as np
class Multiphysics(ph.Physics):
    def getFluxInner(self):
        pass
    def getFluxBoundary(self):
        pass

    def getSourceValue(self, point, volume, variable):
        return 0


class PartialSolution():
    def __init__(self, physics:ph.Physics, polygon:geo.Polygon, varnum) -> None:
        self.physics = physics
        self.polygon = polygon
        self.varnum = varnum
    def getFluxInner(self, point, **kvargs):
        if self.polygon.isPointInside(point):
            return self.physics.getFluxInner(**kvargs)
        else:
            coeff_mid = np.zeros((1,self.varnum))
            coeff_mid[0,kvargs["variable"]] = 1
            coeff_neighbour = np.zeros((1,self.varnum))
            coeff_const = 0 
            return coeff_mid,coeff_neighbour,coeff_const
    def getFluxBoundary(self, point, **kvargs):
        if self.polygon.isPointInside(point):
            return self.physics.getFluxBoundary(**kvargs)
        else:
            coeff_mid = np.zeros((1,self.varnum))
            coeff_mid[0,kvargs["variable"]] = 1
            coeff_neighbour = np.zeros((1,self.varnum))
            coeff_const = 0 
            return coeff_mid,coeff_neighbour,coeff_const

    def getSourceValue(self, point, **kvargs):
        if self.polygon.isPointInside(point):
            return self.physics.getSourceValue(point, **kvargs)
        else:
            return 0


class SequentialCoupling():
    def __init__(self, physics:ph.Physics, polygon1:geo.Polygon, varnum) -> None:
        self.physics = physics
        self.polygon1 = polygon1

        self.varnum = varnum
    def getFluxInner(self, point, parameter, sequence_number, **kvargs):
        if self.polygon.isPointInside(point):
            return self.physics.getFluxInner(**kvargs)
        else:
            coeff_mid = np.zeros((1,self.varnum))
            coeff_mid[0,kvargs["variable"]] = 1
            coeff_neighbour = np.zeros((1,self.varnum))
            coeff_const = 0 
            return coeff_mid,coeff_neighbour,coeff_const
    def getFluxBoundary(self,  point, parameter, sequence_number, **kvargs):
        if self.polygon.isPointInside(point):
            return self.physics.getFluxBoundary(**kvargs)
        else:
            coeff_mid = np.zeros((1,self.varnum))
            coeff_mid[0,kvargs["variable"]] = 1
            coeff_neighbour = np.zeros((1,self.varnum))
            coeff_const = 0 
            return coeff_mid,coeff_neighbour,coeff_const

    def getSourceValue(self,  point, parameter, sequence_number, **kvargs):
        if self.polygon.isPointInside(point):
            return self.physics.getSourceValue(point, **kvargs)
        else:
            return 0
        
class ParametricCoupling():
    def __init__(self, physics:ph.Physics, polygon:geo.Polygon,param, meshgrid,  varnum) -> None:
        self.physics = physics
        self.polygon = polygon
        self.varnum = varnum
        self.param = param
    def getFluxInner(self, point, **kvargs):
        return self.physics.getFluxInner(param = param_here, **kvargs)
    def getFluxBoundary(self, point, **kvargs):
        return self.physics.getFluxBoundary(param = param_here, **kvargs)

    def getSourceValue(self, point, **kvargs):
        return self.physics.getSourceValue(param = param_here, **kvargs)

    def getParamAtPoint(self, point):
        self.param