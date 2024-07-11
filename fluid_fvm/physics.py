import numpy as np
import fluid_fvm.geometry as geo

class Physics():
    def getFluxInner(self):
        pass
    def getFluxBoundary(self):
        pass

    def getSourceValue(self, point, volume, variable):
        return 0

    def _Grad(self, neighbour_vector, num_variables=1):
        # Calculates the normalized flux in the middle of the neighbour_vector (d) (grad(Phi)|d) in the given direction for PDE-s in the form of div(grad(Phi)) = c

        # grad(Phi)|d = Vector(  (phi_e-phi_c)/diff_x, (phi_e-phi_c)/diff_y ) =  Fc*phi_c+Fe*phi_e  
        # (diffx and y are the x and y components of the vector pointing from the current to the neighbouring face)
        # If diff_x or y is 0 the respective component of the gradient is zero
        if neighbour_vector[0,0] == 0:
            comp_x = 0
        else:
            comp_x = 1/neighbour_vector[0,0]
        

        if neighbour_vector[1,0] == 0:
            comp_y = 0
        else:
            comp_y = 1/neighbour_vector[1,0]

        #THIS SHOULD BE num by num, but whatever
        Fe = np.zeros((num_variables, 2))
        Fc = np.zeros((num_variables, 2))

        Fe[:,0] = -comp_x
        Fe[:,1] = -comp_y
        
        Fc[:,0] = comp_x
        Fc[:,1] = comp_y
        return Fe,Fc
    
    
    def _Value(self, neighbour_vector, num_variables = 1):
        # Calculates the normalized x component of the Phi*A_norm vector at the middle of the normal vector, where A_norm is the norml vector of the surface in the given direction.
        # Used for PDE-s in the form of div(Phi) = c

        # Phi|d*(A_norm*e_x) = (phi_c+phi_e)/2*(A_norm*e_x) = Fc*phi_c+Fe*phi_e where e_x is the unit vector in the x direction

        # normalize(A_norm*e_x) = sgn(A_norm.x*e_x.x) = sgn(A_norm.x) (-1, if x is negative, 0 if x is 0 and 1 if x positive)
        mrx = np.eye(num_variables )
        return mrx*1/2, mrx*1/2

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
    

    def getFluxInner(self, material, face_normal,  neighbour_vector, variable):
        gamma = material.getProperty("gamma")

        neighbour_vector_np = neighbour_vector.toNpArray()
        face_normal_np = face_normal.toNpArray()


        Fe, Fc = self._Grad(neighbour_vector=neighbour_vector_np)
        # J = gamma * dot(grad(phi), Se) =gamma*{ (Fc*phi_c+Fe*phi_e).x*Se.x + (Fc*phi_c+Fe*phi_e).y*Se.y } = gamma*{  Jc*phi_c + Je*phi_e } 

        Jc = Fc@face_normal_np
        Je = Fe@face_normal_np
        #print("Jc, Je")
        #print(Jc)
        #print(Je)
        #print("neighbour vetor: "+str(neighbour_vector.x)+ " y "+str(neighbour_vector.y))

        coeff_mid = gamma*Jc
        coeff_neighbour = gamma*Je
        coeff_const = 0

        return coeff_mid, coeff_neighbour, coeff_const
    
    def getFluxBoundary(self, boundary_face_name, material, face_normal,  neighbour_vector, variable):
        gamma = material.getProperty("gamma")

        boundary =self.boundaries[boundary_face_name]
        neighbour_vector_np = neighbour_vector.toNpArray()
        face_normal_np = face_normal.toNpArray()

        if boundary.type =="Dirichlet":
            Fb, Fc = self._Grad(neighbour_vector=neighbour_vector_np)
            # J = gamma * dot(grad(phi), Sb) =gamma*{ (Fc*phi_c+Fb*phi_b).x*Sb.x + (Fc*phi_c+Fb*phi_b).y*Sb.y } = gamma*{  Jc*phi_c + Je*phi_b }

            Jc = Fc@face_normal_np
            Jb = Fb@face_normal_np

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
    
class CouetteFlow(Physics):
    def __init__(self, assembly, boundaries, flowDirectionUnitVector) -> None:
        super().__init__()
        self.flowDirection = flowDirectionUnitVector
        self.flowNormal = geo.Vector(flowDirectionUnitVector.y, -flowDirectionUnitVector.x)
        self.assembly = assembly
        self.defineBoundaries(boundaries)

    def defineBoundaries(self, bcdict):
        assert(set(bcdict.keys()).issubset(self.assembly.getLineNames()))
        assert(set(bcdict.keys()).issuperset(self.assembly.getBaseLineNames()))
        self.boundaries = bcdict
    def isBoundary(self, name):
        return name in self.boundaries.keys()
    

    def getFluxInner(self, material, face_normal,  neighbour_vector, variable):
        mu = material.getProperty("mu")


        neighbour_vector_np = neighbour_vector.toNpArray()
        face_normal_np = face_normal.toNpArray()

     


        Fe, Fc = self._Grad(neighbour_vector=neighbour_vector_np)
        # J = mu * dot(grad(phi), Se) =gamma*{ (Fc*phi_c+Fe*phi_e).x*Se.x + (Fc*phi_c+Fe*phi_e).y*Se.y } = gamma*{  Jc*phi_c + Je*phi_e } 
        face_dir_corr = (face_normal*self.flowNormal)
        
        Jc = Fc@face_normal_np*abs(face_dir_corr)
        Je = Fe@face_normal_np*abs(face_dir_corr)
        #print("Jc, Je")
        #print(Jc)
        #print(Je)
        #print("neighbour vetor: "+str(neighbour_vector.x)+ " y "+str(neighbour_vector.y))

        coeff_mid = mu*Jc
        coeff_neighbour = mu*Je
        coeff_const = 0

        return coeff_mid, coeff_neighbour, coeff_const
    
    def getFluxBoundary(self, boundary_face_name, material, face_normal,  neighbour_vector, variable):
        mu = material.getProperty("mu")

        boundary =self.boundaries[boundary_face_name]
        neighbour_vector_np = neighbour_vector.toNpArray()
        face_normal_np = face_normal.toNpArray()

        if boundary.type =="Dirichlet":
            if abs(face_normal*self.flowNormal) <1e-3:
                raise AttributeError("Dirichlet boundaries must be parallel to flow")
            Fb, Fc = self._Grad(neighbour_vector=neighbour_vector_np)
            # J = gamma * dot(grad(phi), Sb) =gamma*{ (Fc*phi_c+Fb*phi_b).x*Sb.x + (Fc*phi_c+Fb*phi_b).y*Sb.y } = gamma*{  Jc*phi_c + Je*phi_b }

            Jc = Fc@face_normal_np
            Jb = Fb@face_normal_np

            # Division by two due to the fact that the boundary is half as close as the mirrored node on the other side of the boundary
            # Is taken care of in the mesh class
            coeff_mid = mu*Jc
            coeff_neighbour = 0
            coeff_const = -mu*Jb[0]*boundary.value

            return coeff_mid, coeff_neighbour, coeff_const
        elif boundary.type == "Neumann":
            if abs(face_normal*self.flowNormal) >1e-3:
                raise AttributeError("Neumann boundaries must be perpendicular to flow")
            face_length = np.sqrt(face_normal.x**2+face_normal.y**2)
            coeff_mid = 0
            coeff_neighbour = 0
            coeff_const =  boundary.value * face_length
            return coeff_mid, coeff_neighbour, coeff_const
        else:
            raise AttributeError("Boundary condition not supported")
    
class PoissonFlowFixP(Physics):
    def __init__(self, assembly, boundaries, flowDirectionUnitVector, dpdx) -> None:
        super().__init__()
        self.flowDirection = flowDirectionUnitVector
        self.flowNormal = geo.Vector(flowDirectionUnitVector.y, -flowDirectionUnitVector.x)
        self.assembly = assembly
        self.defineBoundaries(boundaries)
        self.dpdx = dpdx

    def defineBoundaries(self, bcdict):
        assert(set(bcdict.keys()).issubset(self.assembly.getLineNames()))
        assert(set(bcdict.keys()).issuperset(self.assembly.getBaseLineNames()))
        self.boundaries = bcdict
    def isBoundary(self, name):
        return name in self.boundaries.keys()
    

    def getFluxInner(self, material, face_normal,  neighbour_vector, variable):
        neighbour_vector_np = neighbour_vector.toNpArray()
        face_normal_np = face_normal.toNpArray()
        mu = material.getProperty("mu")



        Fe, Fc = self._Grad(neighbour_vector=neighbour_vector_np)
        # J = mu * dot(grad(phi), Se) =gamma*{ (Fc*phi_c+Fe*phi_e).x*Se.x + (Fc*phi_c+Fe*phi_e).y*Se.y } = gamma*{  Jc*phi_c + Je*phi_e } 
        face_dir_corr = (face_normal*self.flowNormal)
        
        Jc = Fc@face_normal_np*abs(face_dir_corr)
        Je = Fe@face_normal_np*abs(face_dir_corr)
        #print("Jc, Je")
        #print(Jc)
        #print(Je)
        #print("neighbour vetor: "+str(neighbour_vector.x)+ " y "+str(neighbour_vector.y))

        coeff_mid = mu*Jc
        coeff_neighbour = mu*Je
        coeff_const = 0 #face_dir_corr*self.dpdx*volume

        return coeff_mid, coeff_neighbour, coeff_const
    
    def getFluxBoundary(self, boundary_face_name, material, face_normal,  neighbour_vector, variable):
        mu = material.getProperty("mu")

        boundary =self.boundaries[boundary_face_name]
        neighbour_vector_np = neighbour_vector.toNpArray()
        face_normal_np = face_normal.toNpArray()
        if boundary.type =="Dirichlet":
            if abs(face_normal*self.flowNormal) <1e-3:
                raise AttributeError("Dirichlet boundaries must be parallel to flow")
            Fb, Fc = self._Grad(neighbour_vector=neighbour_vector_np)
            # J = gamma * dot(grad(phi), Sb) =gamma*{ (Fc*phi_c+Fb*phi_b).x*Sb.x + (Fc*phi_c+Fb*phi_b).y*Sb.y } = gamma*{  Jc*phi_c + Je*phi_b }

            Jc = Fc@face_normal_np
            Jb = Fb@face_normal_np

            # Division by two due to the fact that the boundary is half as close as the mirrored node on the other side of the boundary
            # Is taken care of in the mesh class
            coeff_mid = mu*Jc
            coeff_neighbour = 0
            coeff_const = -mu*Jb[0]*boundary.value

            return coeff_mid, coeff_neighbour, coeff_const
        elif boundary.type == "Neumann":
            if abs(face_normal*self.flowNormal) >1e-3:
                raise AttributeError("Neumann boundaries must be perpendicular to flow")
            face_length = np.sqrt(face_normal.x**2+face_normal.y**2)
            coeff_mid = 0
            coeff_neighbour = 0
            coeff_const =  boundary.value * face_length
            return coeff_mid, coeff_neighbour, coeff_const
        else:
            raise AttributeError("Boundary condition not supported")
    
    def getSourceValue(self, point, volume, variable):
        return self.dpdx*volume
    

class Boundary():
    def __init__(self, type, value) -> None:
        self.type = type
        self.value = value
