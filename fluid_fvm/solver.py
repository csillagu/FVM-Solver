import numpy as np
import scipy as sp

class Solver():
    pass


class FullBlockSolver(Solver):
    def __init__(self, Amrx, Bmrx) -> None:
        super().__init__()
        self.Amrx = Amrx
        self.Bmrx = Bmrx

    def solve(self):
        self.res = np.linalg.solve(self.Amrx, self.Bmrx)

class DualSequentialSolver(Solver):
    def __init__(self,sequence,  discretizer1, solver1, discretizer2, solver2) -> None:
        super().__init__()
        self.discretizer1 = discretizer1
        self.sequence = sequence
        self.solver1 = solver1

        self.discretizer2 = discretizer2
        self.solver2 = solver2
    
    def solve(self):
        disc1 = self.discretizer1.discretize(0, self.sequence[0])
        self.solver1.solve(disc1.Amrx, disc1.Bmrx)
        self.res1= self.solver1.res

        disc2 = self.discretizer1.discretize(0, self.sequence[0], self.res1)
        self.solver2.solve(disc2.Amrx, disc2.Bmrx)
        self.res2= self.solver2.res

