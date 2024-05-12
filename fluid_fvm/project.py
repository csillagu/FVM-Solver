
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
    


class Material():
    def __init__(self, name, **kvargs) -> None:
        self.properties = kvargs
        self.name = name

    def getProperty(self, name):
        if name in self.properties:
            return self.properties[name]
        else:
            return None