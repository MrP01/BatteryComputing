# Python code for solving chronoamperometre problem

import numpy as np
from scipy.sparse import diags
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# from scipy.sparse import dia_matrix, dok_matrix
from scipy.sparse.linalg import spsolve

# Python code for External viewer on Mac
plt.switch_backend("MacOSX")
np.set_printoptions(linewidth=200)

# Function creating the scheme matrix for "a"
# Admittedly I used a for loop which was not necessary.
def schemeD(Nx, Nt, CFL, Dirichlet):
    # diagonal vector with (Nx+1)*(Nt+1) entries
    diag = np.ones(Nx * Nt) * (1 + 2 * CFL)

    # off diagonal vectors
    udiag = np.ones(Nx * Nt) * (-CFL)
    ldiag = np.ones(Nx * Nt) * (-CFL)

    # lonely diagonal vector accounts for the explicit coefficient
    lonelydiag = np.ones(Nx * Nt - Nx) *(-1)

    # We aim to zero the entries of the lonely vector corresponding to previous timestep
    # We aim to have tridiagonal matrix (lowerdiag, diag, upperdiag) with (-CFl, 1+2CFL, CFL)
    for t in range(0, Nt):
        first = t * Nx
        last = first + Nx - 1

        diag[first] = 1
        diag[last] = 1

        ldiag[last] = 0  # ultimate
        ldiag[last - 1] = 0  # penultimate

        udiag[first] = 0 if Dirichlet else -1
        udiag[last] = 0

        if t != Nt - 1:
            lonelydiag[first] = 0
            lonelydiag[last] = 0

    # initial condition:
    for k in range(0, Nx):
        diag[k] = 1
        ldiag[k] = 0
        udiag[k] = 0

    ldiag = ldiag[:-1]
    udiag = udiag[:-1]

    # This creates sparse matrix out of diagonals
    return diags([diag, udiag, ldiag, lonelydiag], [0, 1, -1, -Nx], format="csr")

def schemeVolta(Nx, Nt, CFL, isA):
    # diagonal vector with (Nx+1)*(Nt+1) entries
    diag = np.ones(Nx * Nt) * (1 + 2 * CFL)

    # off diagonal vectors
    udiag = np.ones(Nx * Nt) * (-CFL)
    ldiag = np.ones(Nx * Nt) * (-CFL)

    # lonely diagonal vector accounts for the explicit coefficient
    lonelydiag = np.ones(Nx * Nt - Nx) *(-1)

    # second upper diagonal
    u2diag = np.zeros(Nx * Nt)

    # We aim to zero the entries of the lonely vector corresponding to previous timestep
    # We aim to have tridiagonal matrix (lowerdiag, diag, upperdiag) with (-CFl, 1+2CFL, CFL)
    for t in range(0, Nt):
        first = t * Nx
        last = first + Nx - 1

        diag[first] = 1
        diag[last] = 1

        ldiag[last] = 0  # ultimate
        ldiag[last - 1] = 0  # penultimate

        udiag[first] = -1
        udiag[last] = 0

        if t != Nt - 1:
            lonelydiag[first] = 0
            lonelydiag[last] = 0

        if t != 0:
            u2diag[first] = -1/2
            udiag[first] = 2
            diag[first] = -3/2

    # initial condition:
    for k in range(0, Nx):
        diag[k] = 1
        ldiag[k] = 0
        udiag[k] = 0

    ldiag = ldiag[:-1]
    udiag = udiag[:-1]

    # This creates sparse matrix out of diagonals
    return diags([diag, udiag, ldiag, lonelydiag,u2diag], [0, 1, -1, -Nx, 2], format="dia_matrix")


def look_at_time(x_A, time_index, Nx):
    return x_A[time_index * Nx : (time_index + 1) * Nx]


def rhsDirichlet(Nx, Nt, IC, BC1, BC2):
    r = np.zeros(Nx * Nt)
    # Initial conditions for the first Nx rows of r.
    r[0:Nx] = IC

    # Fill the boundary conditions in
    for t in range(1, Nt):
        r[t * Nx] = BC1
        r[t * Nx + Nx - 1] = BC2
    return r


def rhsNeumann(Nx, Nt, IC, BC1, BC2, a, D):
    r = np.zeros(Nx * Nt)
    # Initial conditions for the first Nx rows of r.
    r[0:Nx] = IC

    # Fill the bounary conditions in
    for t in range(1, Nt):
        r[t * Nx] = -D * (3/2* a[t * Nx] -2* a[t * Nx + 1] +1/2*a[t*Nx+2])
        r[t * Nx + Nx-1] = 0
    return r

def rhsNeumannChrono(Nx, Nt, IC, BC1, BC2, dA, D):
    r = np.zeros(Nx * Nt)
    # Initial conditions for the first Nx rows of r.
    r[0:Nx] = IC

    # Fill the bounary conditions in
    for t in range(1, Nt):
        r[t * Nx] = 0.5*()
        r[t * Nx + Nx-1] = 0
    return r


def chronoamperometry():
    ## Set up for solving the heat equation for a
    # Boundaries for space and time
    x0 = 0
    xn = 1  # i replaced a by x_0 and b by xn #should be xn=10
    T0 = 0
    T = 0.2

    # Number of meshpoints and meshsizes
    Nx = 300
    Nt = 400
    dx = (xn - x0) / Nx
    dt = (T - T0) / Nt

    # Courant Friedrichs-Lewy
    CFL = dt / dx**2

    # Values at boundary in space and time
    BC1 = 0
    BC2 = 1
    IC_A = 1
    IC_B =0

    # Set up matrices for solving the Dirichlet scheme for "a"
    A = schemeD(Nx, Nt, CFL, True)
    rhsA = rhsDirichlet(Nx, Nt, IC_A, BC1, BC2)
    np.set_printoptions(linewidth=200)
    x_A = spsolve(A, rhsA)

    # Set up matrices for solving the Neumann scheme for "b"
    B = schemeD(Nx, Nt, CFL, False)
    rhsB = rhsNeumann(Nx, Nt, IC_B, BC1, BC2, x_A, 1)
    x_B = spsolve(B, rhsB)

    # Plotting
    fig, ax = plt.subplots(figsize=(5, 3))
    ax.set(xlim=(x0, xn), ylim=(x_B.min(), x_B.max()))
    x = np.linspace(x0, xn, Nx)
    # t = np.linspace(T0, T, Nt)

    line = ax.plot(x, x_A[0:Nx], color="k", lw=2)[0]  # x_A(0), ...x_A(Nx-1)
    line2 = ax.plot(x, x_B[0:Nx], color="r", lw=2)[0]  # x_A(0), ...x_A(Nx-1)

    def animate(t):
        first = t * Nx
        last = first + Nx - 1
        line.set_ydata(x_A[first : last + 1])
        line2.set_ydata(x_B[first: last + 1])

    anim = FuncAnimation(fig, animate, interval=dt * 4000, frames=Nt )
    anim.save("new.gif", writer="imagemagick")
    plt.show()
    return x_A,x_B

def voltametry():
    ## Set up for solving the heat equation for a
    # Boundaries for space and time
    x0 = 0
    xn = 1  # i replaced a by x_0 and b by xn #should be xn=10
    T0 = 0
    T = 0.2

    # Number of meshpoints and meshsizes
    Nx = 300
    Nt = 400
    dx = (xn - x0) / Nx
    dt = (T - T0) / Nt

    # Courant Friedrichs-Lewy
    CFL = dt / dx**2

    # Values at boundary in space and time
    BC1 = 0
    BC2 = 1
    IC_A = 1
    IC_B =0

    # Set up matrices for solving the Neumann scheme for "a".
    # Input true if a, false if b
    A = schemeVolta(Nx, Nt, CFL, True)
    rhsA = rhsDirichlet(Nx, Nt, IC_A, BC1, BC2)
    np.set_printoptions(linewidth=200)
    x_A = spsolve(A, rhsA)

    # Set up matrices for solving the Neumann scheme for "b"
    B = schemeVolta(Nx, Nt, CFL, False)
    rhsB = rhsNeumann(Nx, Nt, IC_B, BC1, BC2, x_A, 1)
    x_B = spsolve(B, rhsB)

    #Initial values
    E_start = -10
    E_0 = 0
    t_rev=20
    kappa=35
    alpha=1/2

    def Pot(t):
        E_start + t if t<t_rev else E_start - t+2*t_rev

if __name__ == "__main__":
    [x_A,x_B] = chronoamperometry()


