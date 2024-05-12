
class Component():

    def __init__(self) -> None:
        pass

    def setAssembly(self, a):
        self.assembly = a


class Material():
    def __init__(self, **kvargs) -> None:
        self.properties = kvargs
    
    def getProperty(self, name):
        return self.properties[name]