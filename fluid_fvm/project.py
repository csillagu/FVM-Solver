
class Component():

    def __init__(self) -> None:
        self.materials = []
        self.material_connectivity ={}

    def setAssembly(self, a):
        self.assembly = a

    def assignMaterial(self, name, material):
        if material.name not in [mat.name for mat in self.materials]:
            self.materials.append(material)
        self.material_connectivity[name] = material.name
    
    def plot(self, ax):
        colorMap = dict( (p_name,  self._findMaterialByName( self.material_connectivity[p_name]).getProperty("color")) for p_name in self.material_connectivity )
        #Poly1 : mat1.name, poly2_mat2.name
        print(colorMap)
        self.assembly.plot(ax,colorMap=colorMap)

    def _findMaterialByName(self, name):
        for mat in self.materials:
            if mat.name == name:
                return mat
        return None


class Material():
    def __init__(self, name, **kvargs) -> None:
        self.properties = kvargs
        self.name = name
        if self.getProperty("color") is None:
            self.properties["color"] = "#c8ffc8"

    def getProperty(self, name):
        if name in self.properties:
            return self.properties[name]
        else:
            return None