
import numpy as np

import fluid_fvm.mesh as ms
import fluid_fvm.physics as ph
import fluid_fvm.geometry as geo
import fluid_fvm.project as pr

class Discretizer:
    def __init__(self) -> None:
        pass


class LinearFullDiscretizer(Discretizer):
    def __init__(self, component: pr.Component, num_variables) -> None:
        super().__init__()
        
        self.component = component
        self.num_variables = num_variables

        self.element_num = self.component.mesh.getVolumeNodeNum()
        self.Amrx = np.zeros((self.element_num*self.num_variables, self.element_num*self.num_variables))
        self.Bmrx = np.zeros((self.element_num*self.num_variables, 1))
        

    def discretize(self):
        physics: ph.Physics = self.component.physics
        mesh: ms.MeshConfig = self.component.mesh
        geometry: geo.Assembly =self.component.assembly

        for node in range(self.element_num):
            node_point = mesh.getVNode(node)
            node_material = self.component.findMaterialAtPoint(node_point)
            
            neighbour_nodes = mesh.getNeigbouringVolumeVectors(node)
            neighbour_faces = mesh.getNeighbouringFaceLines(node)

            neighbour_node_nums = mesh.getNeighbouringVolumes(node)
            node_volume = mesh.getAreaOfElement(node)
            row_idx = node*self.num_variables
            

            for variable in range(self.num_variables):
                self_coefficient = 0.0
                neighbour_coefficients = np.zeros((len(neighbour_node_nums),self.num_variables), dtype="float64")
                local_B_coefficient = 0.0
                for neighbour in range(len(neighbour_node_nums)):
                    neighbour_line_name = geometry.getCoincidentLineName(neighbour_faces[neighbour])
                    if neighbour_node_nums[neighbour] == []:
                        d_coeffs = physics.getFluxBoundary(material=node_material, face_normal=neighbour_faces[neighbour].getNormal(), 
                                                        neighbour_vector=neighbour_nodes[neighbour], 
                                                        boundary_face_name=neighbour_line_name, variable = variable)
                    else:
                        
                        d_coeffs = physics.getFluxInner(material=node_material, face_normal=neighbour_faces[neighbour].getNormal(), 
                                                        neighbour_vector=neighbour_nodes[neighbour], variable = variable)
                    print("Node:"+str(node) +"  Neighbour: "+str(neighbour) + " Nwighbour idx: " + str(neighbour_node_nums[neighbour]))
                    print("Self"+str(d_coeffs[0]))
                    print("Neighbour:"+str(d_coeffs[1]))
                    print("Const:"+str(d_coeffs[2]))
                    print("---------------")
                    self_coefficient += d_coeffs[0]
                    neighbour_coefficients[neighbour, :] = neighbour_coefficients[neighbour, :]+ d_coeffs[1]
                    local_B_coefficient += d_coeffs[2]

                for neighbour in range(len(neighbour_node_nums)):
                    # Only save the neighbour coefficient into the A matrix for those neighbours, that are not boundaries
                    if neighbour_node_nums[neighbour] != []:
                        neighbour_starting_num =neighbour_node_nums[neighbour]*self.num_variables
                        # print(neighbour_starting_num)
                        # print(neighbour_starting_num+self.num_variables)
                        self.Amrx[row_idx+variable, neighbour_starting_num:(neighbour_starting_num+self.num_variables)] = neighbour_coefficients[neighbour, :]
                        
                #print(neighbour_coefficients)
                #print(str(row_idx)+":"+str(row_idx+self.num_variables)+"<="+str(self_coefficient))
                #print(self.Amrx[row_idx+variable, (row_idx):(row_idx+self.num_variables)])
                self.Amrx[row_idx+variable, (row_idx):(row_idx+self.num_variables)] = self_coefficient
                self.Bmrx[row_idx+variable,0] = local_B_coefficient+physics.getSourceValue(point=node_point, volume = node_volume, variable=variable)
            
