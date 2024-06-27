
import numpy as np

import fluid_fvm.mesh as ms
import fluid_fvm.physics as ph
import fluid_fvm.geometry as geo
import fluid_fvm.project as pr

class Discretizer:
    def __init__(self) -> None:
        pass


class LinearFullDiscretizer(Discretizer):
    def __init__(self, component: pr.Component) -> None:
        super().__init__()
        
        self.component = component
        self.element_num = self.component.mesh.getVolumeNodeNum()
        self.Amrx = np.zeros((self.element_num, self.element_num))
        self.Bmrx = np.zeros((self.element_num, 1))
        

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

            self_coefficient = 0
            neighbour_coefficients = np.zeros((len(neighbour_node_nums),1))
            local_B_coefficient = 0
            
            for neighbour in range(len(neighbour_node_nums)):
                neighbour_line_name = geometry.getCoincidentLineName(neighbour_faces[neighbour])
                if neighbour_node_nums[neighbour] == []:
                    d_coeffs = physics.getFluxBoundary(material=node_material, face_normal=neighbour_faces[neighbour].getNormal(), neighbour_vector=neighbour_nodes[neighbour], boundary_face_name=neighbour_line_name)
                else:
                    #print("Node:"+str(node) +"  Neighbour: "+str(neighbour) + " Nwighbour idx: " + str(neighbour_node_nums[neighbour]))
                    d_coeffs = physics.getFluxInner(material=node_material, face_normal=neighbour_faces[neighbour].getNormal(), neighbour_vector=neighbour_nodes[neighbour])
                self_coefficient += d_coeffs[0]
                neighbour_coefficients[neighbour] += d_coeffs[1]
                local_B_coefficient += d_coeffs[2]

            for neighbour in range(len(neighbour_node_nums)):
                self.Amrx[node, neighbour_node_nums[neighbour]] = neighbour_coefficients[neighbour]

            self.Amrx[node, node] = self_coefficient
            self.Bmrx[node,0] = local_B_coefficient
            
