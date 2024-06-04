import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

def create_sparse_matrix(stencil, n):
    """
    Creates sparse matrix given stencil and length n 

    """
    stencil = np.array(stencil)
    m = len(stencil)
    diagonals = []

    # creeate diagonal for matrix
    for i in range(m):
        diagonals.append(np.full(n, stencil[i]))

    #positioning of the diagonal (for boundaries etc)
    offsets = np.arange(-(m//2), m//2 + 1)
    
    # Create sparse Matrix and safe in "better" format (CSR)
    A = sp.diags(diagonals, offsets, shape=(n, n), format='csr')
    
    return A

def solve_system(stencil, b):
    """
    Solvesm A*x = b based on stencil and Vektor b.

    """
    n = len(b)
    A = create_sparse_matrix(stencil, n)
    x = spla.spsolve(A, b)
    return x


def main():
    stencil = [1, -2, 1]
    b = np.array([1, 2, 3, 4, 5])

    x = solve_system(stencil, b)
    print(x)
    pass

if __name__ == "__main__":
    main()