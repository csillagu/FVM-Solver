
import numpy as np


class Discretizer:
    def __init__(self) -> None:
        pass


class LinearFullDiscretizer(Discretizer):
    def __init__(self, element_size) -> None:
        super().__init__()
        self.element_size = element_size
        self.Amrx = np.zeros((element_size, element_size))
        self.Bmrx = np.zeros((element_size, 1))

    def discretize(self, physics, mesh, materials, geometry):
        pass
