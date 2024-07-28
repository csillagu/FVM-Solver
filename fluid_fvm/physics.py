import numpy as np
import fluid_fvm.geometry as geo

class Physics():
    def __init__(self) -> None:
        self.e_x = np.array([1,0])
        self.e_y = np.array([[0, 1]])
        
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


class HeatTransferFull(Physics):
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
    

    def getFluxInner(self, param, material, face_normal,  neighbour_vector, variable, ):
        gamma = material.getProperty("gamma")
        velocity = np.transpose(param.toNpArray())

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

        # (grad T)[x] = sum(T|f * face_normal*1/vol*e_x)
        Fe_gradT, Fc_gradT = self._Value(neighbour_vector=neighbour_vector_np)

        Jc_gradT = Fc_gradT*velocity@face_normal_np*1/1
        Je_gradT = Fe_gradT*velocity@face_normal_np*1/1


        coeff_mid = gamma*Jc+Jc_gradT
        coeff_neighbour = gamma*Je+Je_gradT
        coeff_const = 0

        return coeff_mid, coeff_neighbour, coeff_const
    
    def getFluxBoundary(self, param,boundary_face_name, material, face_normal,  neighbour_vector, variable, ):
        velocity = np.transpose(param.toNpArray())
        gamma = material.getProperty("gamma")

        boundary =self.boundaries[boundary_face_name]
        neighbour_vector_np = neighbour_vector.toNpArray()
        face_normal_np = face_normal.toNpArray()

        if boundary.type =="Dirichlet":
            # T = const
            Fb, Fc = self._Grad(neighbour_vector=neighbour_vector_np)
            # J = gamma * dot(grad(phi), Sb) =gamma*{ (Fc*phi_c+Fb*phi_b).x*Sb.x + (Fc*phi_c+Fb*phi_b).y*Sb.y } = gamma*{  Jc*phi_c + Je*phi_b }

            Jc = Fc@face_normal_np
            Jb = Fb@face_normal_np

            # (grad p)[x] = sum(p|f * face_normal*1/vol*e_x), here p|f = boundary.value

            Jconst = boundary.value*velocity@face_normal_np*1/1

            # Division by two due to the fact that the boundary is half as close as the mirrored node on the other side of the boundary
            # Is taken care of in the mesh class
            coeff_mid = gamma*Jc
            coeff_neighbour = 0
            coeff_const = -gamma*Jb*boundary.value-Jconst

            return coeff_mid, coeff_neighbour, coeff_const
        elif boundary.type == "Neumann":
            face_length = np.sqrt(face_normal.x**2+face_normal.y**2)

            # (grad T)[x] = sum(T|f * face_normal*1/vol*velocity)
            Jc_gradT = 1*velocity@face_normal_np*1/1 #T|f = T1 thus coefficient 1

            # Division by two due to the fact that the boundary is half as close as the mirrored node on the other side of the boundary
            # Is taken care of in the mesh class
            coeff_mid = Jc_gradT
            coeff_neighbour = 0
            coeff_const =  boundary.value * face_length
            return coeff_mid, coeff_neighbour, coeff_const
        else:
            raise AttributeError("Boundary condition not supported")
    def getSourceValue(self, param,point, volume, variable, ):
        return 0

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

class StokesFlow(Physics):
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
        if variable ==0:
            return self.getFluxInnerMomentumX( material, face_normal,  neighbour_vector)
        if variable == 1:
            return self.getFluxInnerMomentumY( material, face_normal,  neighbour_vector)
        if variable == 2:
            return self.getFluxInnerMassConserv( material, face_normal,  neighbour_vector)
        
    def getFluxInnerMomentumX(self, material, face_normal,  neighbour_vector):
        neighbour_vector_np = neighbour_vector.toNpArray()
        face_normal_np = face_normal.toNpArray()
    
        # 0 = -div p +mu*div grad v
        mu = material.getProperty("mu")


        # (div grad u)[x] = sum (grad u |f * face_normal * 1/Vol *e_x)
        Fe_ux, Fc_ux = self._Grad(neighbour_vector=neighbour_vector_np, num_variables=2)
        Je_ux = self.e_x@Fe_ux@face_normal_np*1/1
        Jc_ux = self.e_x@Fc_ux@face_normal_np*1/1
        # (grad p)[x] = sum(p|f * face_normal*1/vol*e_x)
        Fe_p, Fc_p = self._Value(neighbour_vector=neighbour_vector_np)

        Jc_p = Fc_p*self.e_x@face_normal_np*1/1
        Je_p = Fe_p*self.e_x@face_normal_np*1/1

        # mu*div grad u - grad p = 0
        coeff_mid = np.array([mu*Jc_ux[0], 0,Jc_p[0,0]])
        coeff_neighbour = np.array([mu*Je_ux[0],0, Je_p[0,0]])
        coeff_const = 0 

        return coeff_mid, coeff_neighbour, coeff_const
    
    def getFluxInnerMomentumY(self, material, face_normal,  neighbour_vector):
        # 0 = -div p +mu*div grad v
        neighbour_vector_np = neighbour_vector.toNpArray()
        face_normal_np = face_normal.toNpArray()
        mu = material.getProperty("mu")


        # (div grad u)[y] = sum (grad u |f * face_normal * 1/Vol *e_y)
        Fe_uy, Fc_uy = self._Grad(neighbour_vector=neighbour_vector_np, num_variables=2)
        Je_uy = self.e_y@Fe_uy@face_normal_np*1/1
        Jc_uy = self.e_y@Fc_uy@face_normal_np*1/1
        # (grad p)[y] = sum(p|f * face_normal*1/vol*e_y)
        Fe_p, Fc_p = self._Value(neighbour_vector=neighbour_vector_np)

        Jc_p = Fc_p*self.e_y@face_normal_np*1/1
        Je_p = Fe_p*self.e_y@face_normal_np*1/1

        # mu*div grad u - grad p = 0
        # coeff_u, coeff_v, coeff_p
        coeff_mid = np.array([0, mu*Jc_uy[0,0], Jc_p[0,0]])
        coeff_neighbour = np.array([0, mu*Je_uy[0,0], Je_p[0,0]])
        coeff_const = 0 

        return coeff_mid, coeff_neighbour, coeff_const
    
    def getFluxInnerMassConserv(self, material, face_normal, neighbour_vector):
        # div v = 0
        neighbour_vector_np = neighbour_vector.toNpArray()
        face_normal_np = face_normal.toNpArray()
        # div v = dv.x/dx+dv.y/dy = sum(v|f * face_normal*1/vol)
        Fe_u, Fc_u = self._Value(neighbour_vector=neighbour_vector_np, num_variables = 2)
        Jc_u = Fc_u@face_normal_np*1/1
        Je_u = Fe_u@face_normal_np*1/1

        #print("Jc, Je")
        #print(Jc)
        #print(Je)
        #print("neighbour vetor: "+str(neighbour_vector.x)+ " y "+str(neighbour_vector.y))

        
        coeff_mid = np.array([Jc_u[0,0], Jc_u[1,0],0])
        coeff_neighbour = np.array([Je_u[0,0], Je_u[1,0],0])
        coeff_const = 0 

        return coeff_mid, coeff_neighbour, coeff_const
    
    def getFluxBoundary(self, boundary_face_name, material, face_normal,  neighbour_vector, variable):
        if variable ==0:
            return self.getFluxBoundaryMomentumX( boundary_face_name, material, face_normal,  neighbour_vector)
        if variable ==1:
            return self.getFluxBoundaryMomentumY( boundary_face_name, material, face_normal,  neighbour_vector)
        if variable == 2:
            return self.getFluxBoundaryMassConserv( boundary_face_name, material, face_normal,  neighbour_vector)

    def getFluxBoundaryMomentumX(self, boundary_face_name, material, face_normal,  neighbour_vector):
        mu = material.getProperty("mu")
        neighbour_vector_np = neighbour_vector.toNpArray()
        face_normal_np = face_normal.toNpArray()
        boundary =self.boundaries[boundary_face_name]

        if boundary.type =="Velocity_inlet":

            Fb_u, Fc_u = self._Grad(neighbour_vector=neighbour_vector_np, num_variables=2)
            # (div grad u)[x] = sum (grad u |f * face_normal * 1/Vol *e_x)

            Jc_u = self.e_x@Fc_u@face_normal_np*1/1
            Jb_u = self.e_x@Fb_u@face_normal_np*1/1

            # (grad p)[x] = sum(p|f * face_normal*1/vol*e_x)
            Jc_p = 1*self.e_x@face_normal_np*1/1

            # Division by two due to the fact that the boundary is half as close as the mirrored node on the other side of the boundary
            # Is taken care of in the mesh class
            coeff_mid = np.array([mu*Jc_u[0], 0, Jc_p[0]])
            coeff_neighbour = np.array([0,0, 0])
            coeff_const = -mu*Jb_u[0]*boundary.value # maybe -??

            return coeff_mid, coeff_neighbour, coeff_const
        elif boundary.type == "Pressure_outlet":
            # Directional derivative is 0, thus grad v|f*face_normal = 0, to achieve this we choose grad v|f = 0
            Fb_u = np.zeros((2,2))
            Fc_u = np.zeros((2,2))

            # (div grad u)[x] = sum (grad u |f * face_normal * 1/Vol *e_x)
            Jc_u = self.e_x@Fc_u@face_normal_np*1/1
            Jb_u = self.e_x@Fb_u@face_normal_np*1/1

            # (grad p)[x] = sum(p|f * face_normal*1/vol*e_x), here p|f = boundary.value

            Jconst = boundary.value*self.e_x@face_normal_np*1/1

            coeff_mid = np.array([mu*Jc_u[0], 0,0])
            coeff_neighbour = np.array([0,0, 0])
            coeff_const = Jconst[0] 
            return coeff_mid, coeff_neighbour, coeff_const

        elif boundary.type == "No_slip":

            Fb_u, Fc_u = self._Grad(neighbour_vector=neighbour_vector_np, num_variables=2)

            # (div grad u)[x] = sum (grad u |f * face_normal * 1/Vol *e_x)
            Jc_u = self.e_x@Fc_u@face_normal_np*1*1
            Jb_u = self.e_x@Fb_u@face_normal_np*1/1

            # (grad p)[x] = sum(p|f * face_normal*1/vol*e_x)
            Jc_p = 1*self.e_x@face_normal_np*1/1 #p|f = p1 thus coefficient 1

            # Division by two due to the fact that the boundary is half as close as the mirrored node on the other side of the boundary
            # Is taken care of in the mesh class
            coeff_mid = np.array([mu*Jc_u[0], 0, Jc_p[0]])
            coeff_neighbour = np.array([0,0, 0])
            coeff_const = mu*Jb_u[0]*0 # boundary value is 0, as there is no slip

            return coeff_mid, coeff_neighbour, coeff_const
        else:
            raise AttributeError("Boundary condition not supported")
        

    def getFluxBoundaryMomentumY(self, boundary_face_name, material, face_normal,  neighbour_vector):
        mu = material.getProperty("mu")
        neighbour_vector_np = neighbour_vector.toNpArray()
        face_normal_np = face_normal.toNpArray()
        boundary =self.boundaries[boundary_face_name]

        if boundary.type =="Velocity_inlet":

            Fb_u, Fc_u = self._Grad(neighbour_vector=neighbour_vector_np, num_variables=2)
            # (div grad u)[y] = sum (grad u |f * face_normal * 1/Vol *e_y)

            Jc_u = self.e_y@Fc_u@face_normal_np*1/1
            Jb_u = self.e_y@Fb_u@face_normal_np*1/1

            # (grad p)[y] = sum(p|f * face_normal*1/vol*e_y)
            Jc_p = 1*self.e_y@face_normal_np*1/1

            # Division by two due to the fact that the boundary is half as close as the mirrored node on the other side of the boundary
            # Is taken care of in the mesh class
            coeff_mid = np.array([0,mu*Jc_u[0,0], Jc_p[0,0]])
            coeff_neighbour = np.array([0,0, 0])
            coeff_const = mu*Jb_u[0]*0 # V_input is 0

            return coeff_mid, coeff_neighbour, coeff_const
        elif boundary.type == "Pressure_outlet":
            # Directional derivative is 0, thus grad v|f*face_normal = 0, to achieve this we choose grad v|f = 0
            Fb_u = np.zeros((2,2))
            Fc_u = np.zeros((2,2))

            # (div grad u)[y] = sum (grad u |f * face_normal * 1/Vol *e_y)
            Jc_u = Fc_u*self.e_y@face_normal_np*1/1
            Jb_u = Fb_u*self.e_y@face_normal_np*1/1

            # (grad p)[y] = sum(p|f * face_normal*1/vol*e_y), here p|f = boundary.value

            Jconst = boundary.value*self.e_y@face_normal_np*1/1

            coeff_mid = np.array([0, mu*Jc_u[0,0], 0])
            coeff_neighbour = np.array([0,0, 0])
            coeff_const = Jconst[0]  
            return coeff_mid, coeff_neighbour, coeff_const

        elif boundary.type == "No_slip":

            Fb_u, Fc_u = self._Grad(neighbour_vector=neighbour_vector_np, num_variables=2)

            # (div grad u)[y] = sum (grad u |f * face_normal * 1/Vol *e_y)
            Jc_u = self.e_y@Fc_u@face_normal_np*1/1
            Jb_u = self.e_y@Fb_u@face_normal_np*1/1

            # (grad p)[y] = sum(p|f * face_normal*1/vol*e_y)
            Jc_p = self.e_y@face_normal_np*1/1 #p|f = p1 thus coefficient 1

            # Division by two due to the fact that the boundary is half as close as the mirrored node on the other side of the boundary
            # Is taken care of in the mesh class
            coeff_mid = np.array([0, mu*Jc_u[0,0], Jc_p[0,0]])
            coeff_neighbour = np.array([0,0, 0])
            coeff_const = mu*Jb_u[0]*0 # boundary value is 0, as there is no slip

            return coeff_mid, coeff_neighbour, coeff_const
        else:
            raise AttributeError("Boundary condition not supported")
    
    def getFluxBoundaryMassConserv(self, boundary_face_name, material, face_normal,  neighbour_vector):
        mu = material.getProperty("mu")
        neighbour_vector_np = neighbour_vector.toNpArray()
        face_normal_np = face_normal.toNpArray()
        boundary =self.boundaries[boundary_face_name]

        if boundary.type =="Velocity_inlet":

            Fb_u = np.array([boundary.value, 0])
            # div v = dv.x/dx+dv.y/dy = sum(v|f * face_normal*1/vol)
            Jb_u = Fb_u@face_normal_np*1/1
            

            # Division by two due to the fact that the boundary is half as close as the mirrored node on the other side of the boundary
            # Is taken care of in the mesh class
            coeff_mid = np.array([0, 0, 0])
            coeff_neighbour = np.array([0,0, 0])
            coeff_const = -mu*Jb_u[0]

            return coeff_mid, coeff_neighbour, coeff_const
        elif boundary.type == "Pressure_outlet":
            Fb_u = np.eye(2)
            # div v = dv.x/dx+dv.y/dy = sum(v|f * face_normal*1/vol)
            Jb_u = Fb_u@face_normal_np*1/1
            

            # Division by two due to the fact that the boundary is half as close as the mirrored node on the other side of the boundary
            # Is taken care of in the mesh class
            coeff_mid = np.array([mu*Jb_u[0,0],mu*Jb_u[1,0], 0])
            coeff_neighbour = np.array([0,0, 0])
            coeff_const =  0
            return coeff_mid, coeff_neighbour, coeff_const

        elif boundary.type == "No_slip":

            Fb_u = np.array([0, 0])
            # div v = dv.x/dx+dv.y/dy = sum(v|f * face_normal*1/vol)
            Jb_u = Fb_u@face_normal_np*1/1
            

            # Division by two due to the fact that the boundary is half as close as the mirrored node on the other side of the boundary
            # Is taken care of in the mesh class
            coeff_mid = np.array([0, 0, 0])
            coeff_neighbour = np.array([0,0, 0])
            coeff_const = -mu*Jb_u[0]

            return coeff_mid, coeff_neighbour, coeff_const
        else:
            raise AttributeError("Boundary condition not supported")
      
    
class Boundary():
    def __init__(self, type, value) -> None:
        self.type = type
        self.value = value
