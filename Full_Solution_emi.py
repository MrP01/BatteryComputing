# Python code for solving chronoamperometre problem

import numpy as np
from scipy.sparse import diags

# from scipy.sparse import dia_matrix, dok_matrix
from scipy.sparse.linalg import spsolve

# import matplotlib.pyplot as plt
# from matplotlib import cm


# Function creating the scheme matrix for "a"
# Admittedly I used a for loop which was not necessary.
def schemeD(Nx, Nt, CFL, Dirichlet):
    # diagonal vector with (Nx+1)*(Nt+1) entries
    diag = np.ones((Nx + 1) * (Nt + 1)) * (1 + 2 * CFL)

    # off diagonal vectors
    udiag = np.ones((Nx + 1) * (Nt + 1) - 1) * (-CFL)
    ldiag = np.ones((Nx + 1) * (Nt + 1) - 1) * (-CFL)

    # lonely diagonal vector accounts for the explicit coefficient
    lonelydiag = np.ones((Nx + 1) * (Nt + 1) - (Nx + 1))

    # We aim to zero the entries of the lonely vector corresponding to previous timestep
    # We aim to have tridiagonal matrix (lowerdiag, diag, upperdiag) with (-CFl, 1+2CFL, CFL)

    for t in range(0, Nt + 1):
        diag[t * (Nx + 1)] = 1
        diag[t * (Nx + 1) + Nx] = 1

        ldiag[t * (Nx + 1) + Nx - 1] = 0

        if Dirichlet:
            udiag[t * (Nx + 1)] = 0
        else:
            udiag[t * (Nx + 1)] = -1

        if t != Nt:
            ldiag[t * (Nx + 1) + Nx] = 0
            if Dirichlet:
                ldiag[t * (Nx + 1) + Nx] = 0  # to check
            else:
                ldiag[t * (Nx + 1) + Nx - 1] = -1
            udiag[t * (Nx + 1) + Nx] = 0
            lonelydiag[t * (Nx + 1)] = 0
            lonelydiag[t * (Nx + 1) + Nx] = 0

    for k in range(0, Nx + 1):
        diag[k] = 1
        ldiag[k] = 0
        udiag[k] = 0

    # This creates sparse matrix out of diagonals
    A = diags([diag, udiag, ldiag, lonelydiag], [0, 1, -1, -(Nx + 1)],format = 'csr')
    return A


def rhsDirichlet(Nx, Nt, IC, BC1, BC2):
    r = np.zeros((Nx + 1) * (Nt + 1))
    # Initial conditions for the first Nx rows of r.
    r[0:Nx] = IC

    # Fill the bounary conditions in
    for t in range(1, Nt + 1):
        r[t * (Nx + 1)] = BC1
        r[t * (Nx + 1) + Nx] = BC2
    return r


def rhsNeumann(Nx, Nt, IC, BC1, BC2, a):
    r = np.zeros((Nx + 1) * (Nt + 1))
    # Initial conditions for the first Nx rows of r.
    r[0:Nx] = IC

    # Fill the bounary conditions in
    for t in range(1, Nt + 1):
        r[t * (Nx + 1)] = a[0, t] - a[1, t]
        r[t * (Nx + 1) + Nx] = a[-1, t] - a[-2, t]
    return r


if __name__ == "__main__":
    ## Set up for solving the heat equation for a
    # Boundaries for space and time
    x0 = 0
    xn = 10  # i replaced a by x_0 and b by xn #should be xn=10
    T0 = 0
    T = 1

    # Space mesh
    # dx = 0.1 ; dt = 0.01

    # Number of meshpoints and meshsizes
    Nx = 2
    Nt = 2
    dx = (xn - x0) / Nx
    dt = (T - T0) / Nt

    # Courant Friedrichs-Lewy
    CFL = dt / dx**2

    # Values at boundary in space and time
    BC1 = 0
    BC2 = 1
    IC = 1

    # Set up matrices for solving the Dirichlet scheme for "a"
    A = schemeD(Nx, Nt, CFL, True)
    rhsA = rhsDirichlet(Nx, Nt, IC, BC1, BC2)
    print(A.toarray())
    x= spsolve(A,rhsA)
    print("x")
    print(x)



    # Set up matrices for solving the Neumann scheme for "b"
    #B = schemeD(Nx, Nt, CFL, False)
    #rhsB = rhsNeumann(Nx, Nt, IC, BC1, BC2)
