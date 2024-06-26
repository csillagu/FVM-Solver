import numpy as np
import fluid_fvm.geometry as geo

class Physics():
    def getFluxInner(self):
        pass
    def getFluxBoundary(self):
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
        #print("Jc, Je")
        #print(Jc)
        #print(Je)
        #print("neighbour vetor: "+str(neighbour_vector.x)+ " y "+str(neighbour_vector.y))

        coeff_mid = gamma*Jc
        coeff_neighbour = gamma*Je
        coeff_const = 0

        return coeff_mid, coeff_neighbour, coeff_const
    
    def getFluxBoundary(self, boundary_face_name, material, face_normal,  neighbour_vector):
        gamma = material.getProperty("gamma")

        boundary =self.boundaries[boundary_face_name]

        if boundary.type =="Dirichlet":
            Fb, Fc = self._Gradient(neighbour_vector=neighbour_vector)
            # J = gamma * dot(grad(phi), Sb) =gamma*{ (Fc*phi_c+Fb*phi_b).x*Sb.x + (Fc*phi_c+Fb*phi_b).y*Sb.y } = gamma*{  Jc*phi_c + Je*phi_b }

            Jc = Fc*face_normal
            Jb = Fb*face_normal

            # Division by two due to the fact that the boundary is half as close as the mirrored node on the other side of the boundary
            # Is taken care of in the mesh class
            coeff_mid = gamma*Jc
            coeff_neighbour = 0
            coeff_const = -gamma*Jb*boundary.value

            return coeff_mid, coeff_neighbour, coeff_const
        elif boundary.type == "Neumann":
            face_length = np.sqrt(face_normal.x**2+face_normal.y**2)
            coeff_mid = 0
            coeff_neighbour = 0
            coeff_const =  boundary.value * face_length
            return coeff_mid, coeff_neighbour, coeff_const
        else:
            raise AttributeError("Boundary condition not supported")
    
    def _Gradient(self, neighbour_vector):
        # grad(phi) = Vector(  (phi_e-phi_c)/diff_x, (phi_e-phi_c)/diff_y ) =  Fc*phi_c+Fe*phi_e  
        # (diffx and y are the x and y components of the vector pointing from the current to the neighbouring face)
        # If diff_x or y is 0 the respective component of the gradient is zero
        if neighbour_vector.x == 0:
            comp_x = 0
        else:
            comp_x = 1/neighbour_vector.x
        

        if neighbour_vector.y == 0:
            comp_y = 0
        else:
            comp_y = 1/neighbour_vector.y

        Fc = geo.Vector(-(comp_x), -(comp_y))
        Fe = geo.Vector((comp_x), (comp_y))
        return Fc,Fe

class Boundary():
    def __init__(self, type, value) -> None:
        self.type = type
        self.value = value
