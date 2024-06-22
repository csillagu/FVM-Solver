import numpy as np
import fluid_fvm.geometry as geo

class Physics():
    pass

class HeatTransfer(Physics):
    def __init__(self, assembly, boundaries) -> None:
        super().__init__()
        self.assembly = assembly
        self.defineBoundaries(boundaries)

    def defineBoundaries(self, bcdict):
        assert(set(bcdict.keys()).issubset(self.assembly.getLineNames()))
        assert(set(bcdict.keys()).issuperset(self.assembly.getBaseLineNames()))
        self.boundaries = bcdict
    def isBoundary(self, name):
        return name in self.boundaries.keys()
    

    def getFluxInner(self, material, face_normal,  neighbour_vector):
        gamma = material.getProperty("gamma")



        Fe, Fc = self._Gradient(neighbour_vector=neighbour_vector)
        # J = gamma * dot(grad(phi), Se) =gamma*{ (Fc*phi_c+Fe*phi_e).x*Se.x + (Fc*phi_c+Fe*phi_e).y*Se.y } = gamma*{  Jc*phi_c + Je*phi_e } 

        Jc = Fc*face_normal
        Je = Fe*face_normal

        coeff_mid = gamma*Jc
        coeff_neighbour = gamma*Je
        coeff_const = 0

        return coeff_mid, coeff_neighbour, coeff_const
    
    def getFluxBoundary(self, boundary, material, face_normal,  neighbour_vector):
        gamma = material.getProperty("gamma")

        if boundary.type =="Dirichlet":
            Fb, Fc = self._Gradient(neighbour_vector=neighbour_vector)
            # J = gamma * dot(grad(phi), Sb) =gamma*{ (Fc*phi_c+Fb*phi_b).x*Sb.x + (Fc*phi_c+Fb*phi_b).y*Sb.y } = gamma*{  Jc*phi_c + Je*phi_b }

            Jc = Fc*face_normal
            Jb = Fb*face_normal

            # Division by two is due to the fact that the boundary is half as close as the mirrored node on the other side of the boundary
            coeff_mid = gamma*Jc/2 
            coeff_neighbour = 0
            coeff_const = gamma*Jb/2*boundary.value

            return coeff_mid, coeff_neighbour, coeff_const
        else:
            raise AttributeError("Boundary condition not supported")
    
    def _Gradient(self, neighbour_vector):
        # grad(phi) = Vector(  (phi_e-phi_c)/diff_x, (phi_e-phi_c)/diff_y ) =  Fc*phi_c+Fe*phi_e  
        # (diffx and y are the x and y components of the vector pointing from the current to the neighbouring face)
        # If diff_x or y is 0 the respective component of the gradient is zero
        if neighbour_vector.x == 0:
            diff_x = 0.1
            corr_x = 0
        else:
            diff_x = neighbour_vector.x
            corr_x = 1

        if neighbour_vector.y == 0:
            diff_y = 0.1
            corr_y = 0
        else:
            diff_y = neighbour_vector.y
            corr_y = 1  

        Fc = geo.Vector(-(1/diff_x*corr_x), -(1/diff_y*corr_y))
        Fe = geo.Vector((1/diff_x*corr_x), (1/diff_y*corr_y))
        return Fc,Fe

class Boundary():
    def __init__(self, type, value) -> None:
        self.type = type
        self.value = value
