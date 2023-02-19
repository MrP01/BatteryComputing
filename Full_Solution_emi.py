# Python code for solving chronoamperometre problem

import math
import pathlib

from scipy import special
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation
from scipy.sparse import diags
from scipy.sparse.linalg import spsolve

# Python code for External viewer on Mac
#plt.switch_backend("MacOSX")
RESULTS_FOLDER = pathlib.Path(__file__).resolve().parent / "results"

## This function relates to CHRONOAMPEROMETRY  where we first tried to solve for x AND t altogether
# This function creates the scheme MATRIX for "a"
def schemeD(Nx, Nt, muD, Dirichlet):
    # diagonal vector with (Nx+1)*(Nt+1) entries
    diag = np.ones(Nx * Nt) * (1 + 2 * muD)

    # off diagonal vectors
    udiag = np.ones(Nx * Nt) * (-muD)
    ldiag = np.ones(Nx * Nt) * (-muD)

    # lonely diagonal vector accounts for the explicit coefficient
    lonelydiag = np.ones(Nx * Nt - Nx) * (-1)

    # At each timestep, we put the coefficients of the boundary values and their neighbours respectively
    # to 1 and 0. If Neumann, the neighbour coefficient is not 0 but -1
    for t in range(0, Nt):
        first = t * Nx          # first boundary entry of a timestep
        last = first + Nx - 1   # last boundary entry of a timestep

        diag[first] = 1 # first BC entry of a timestep
        diag[last] = 1 # last BC entry of a timestep

        ldiag[last] = 0      # last boundary entry
        ldiag[last - 1] = 0  # penultimate entry

        udiag[first] = 0 if Dirichlet else -1
        udiag[last] = 0

        if t != Nt - 1:
            lonelydiag[first] = 0
            lonelydiag[last] = 0

    # account for initial condition:
    for k in range(0, Nx):
        diag[k] = 1
        ldiag[k] = 0
        udiag[k] = 0

    # shorten the lower diagonal and the upper diagonal to have Nx-1 entries
    ldiag = ldiag[:-1]
    udiag = udiag[:-1]

    # This creates sparse matrix out of specified diagonals
    return diags([diag, udiag, ldiag, lonelydiag], [0, 1, -1, -Nx], format="csr")

# rhsDirichlet assembles the right hand side of matrix system, where we have Dirichlet BC
def rhsDirichlet(Nx, Nt, IC, BC1, BC2):
    r = np.zeros(Nx * Nt)
    # Initial conditions for the first Nx rows of r.
    r[0:Nx] = IC

    # Fill the boundary conditions in
    for t in range(1, Nt):
        r[t * Nx] = BC1
        r[t * Nx + Nx - 1] = BC2
    return r

# rhsNeumann assembles the right hand side of matrix system, where we have Neumann BC
def rhsNeumann(Nx, Nt, IC, BC2, a, D): #Neumann for the first boundary and Dirichlet with the second
    r = np.zeros(Nx * Nt)
    # Initial conditions for the first Nx rows of r.
    r[0:Nx] = IC

    # Fill the boundary conditions in
    for t in range(1, Nt):
        r[t * Nx] = -D * (3 / 2 * a[t * Nx] - 2 * a[t * Nx + 1] + 1 / 2 * a[t * Nx + 2]) # using the Taylor
        r[t * Nx + Nx - 1] = BC2
    return r

## CHRONOAMPEROMETRY is the main code for running experiment with Dirichlet BC for A and Neumann for B
def chronoamperometry():
    ## Set up for solving the heat equation for a
    # Boundaries for space and time
    x0 = 0
    xn = 1  # i replaced a by x_0 and b by xn #should be xn=10
    t0 = 0
    tn = 0.2

    # Number of meshpoints and meshsizes
    Nx = 300
    Nt = 400
    dx = (xn - x0) / Nx
    dt = (tn - t0) / Nt

    # Courant Friedrichs-Lewy
    CFL = dt / dx**2
    D = 1

    # Values at boundary in space and time
    BC1A = 0
    BC2A = 1
    IC_A = 1
    IC_B = 0

    # Set up matrices for solving the Dirichlet scheme for "a"
    A = schemeD(Nx, Nt, CFL * D, True)
    rhsA = rhsDirichlet(Nx, Nt, IC_A, BC1A, BC2A)
    np.set_printoptions(linewidth=200)
    x_A = spsolve(A, rhsA)

    # Set up matrices for solving the Neumann scheme for "b"
    B = schemeD(Nx, Nt, CFL * D, False)
    rhsB = rhsNeumann(Nx, Nt, IC_B, 0, x_A, 1)
    x_B = spsolve(B, rhsB)

    # Store the real solution for a
    real_A = np.zeros(Nx*Nt)
    real_A[0:Nx] = IC_A
    for j in range(1,Nt):
        x=np.linspace(0,xn,Nx)
        h=math.sqrt(j*dt)
        real_A[j*Nx:(j+1)*Nx] = special.erf(x/(2*h))
    print(real_A)

    # Compare real to numerical solutionn
    #plotAnimate(x0,xn,Nx,Nt,dt,x_A,real_A,"chrono_real_sol.gif","Numerical Solution","Real Solution")

    # Plotting the concentration of A versus that of B
    plotAnimate(x0, xn, Nx, Nt, dt, x_A, x_B, "chrono.gif","Concentration A","Concentration B")


# VOLTAMETRY is the main code for running our second experiment, namely using the Butler-Volmer equation
def voltametry():
    ## Set up for solving the heat equation for a
    # Boundaries for space and time
    x0 = 0
    xn = 1
    t0 = 0
    tn = 1

    # Number of meshpoints and meshsizes
    Nx = 300
    dx = (xn - x0) / Nx
    Nt = 400
    dt = (tn - t0) / Nt

    # Courant Friedrichs-Lewy and diffusion constant
    CFL = dt / dx**2
    D = 1

    # Values at boundary in space and time
    IC_A = 1
    IC_B = 0

    # Values found for the Taylor expansion of ∂a/dx at x=0
    alpha = -3 / 2
    beta = 2
    gamma = -1 / 2

    # Initial values
    E_start = -10
    E_0 = 1
    t_rev = 8
    kappa = 35
    alpha = 1 / 2

    def Pot(t):  # This function returns our potential
        E = E_start + t if t < t_rev else E_start - t + 2 * t_rev
        return E

    def G_A(t): # This function represents the first term of our Butler-Volmer function
        return kappa * dx * math.exp((1 - alpha) * (Pot(t) - E_0))

    def G_B(t): # This function represents the positive part of the second term of our Butler-Volmer function
        return kappa * dx * math.exp((-alpha) * (Pot(t) - E_0))

    # oneTimeStepMatrix creates a matrix that incorporates the Butler-Volmer scheme but only for one step.
    # This is similar to schemeD, only that schemeD creates a matrix involving all unknowns in time and space
    def oneTimeStepMatrix(t):
        ## Create a matrix for one timestep. The unknowns are (A0...A_N-1,B0,...B_N-1)
        # The diagonal is 1+2*D*mu except at the boundaries
        diag = np.ones(2 * Nx) * (1 + 2 * D * CFL)
        diag[0] =  G_A(t) +1
        diag[Nx - 1] = 1
        diag[Nx] = -1
        diag[-1] = 1

        # The upper diag is -mu*D except at boundaries
        udiag = np.ones(2 * Nx) * (-D * CFL)
        udiag[0] = -1
        udiag[Nx - 1] = 0
        udiag[Nx] = 1
        udiag = udiag[:-1]  # discard the -1 entry

        # The lower diag is -mu*D except at boundaries
        ldiag = np.ones(2 * Nx) * (-D * CFL)
        ldiag[Nx - 1] = 0
        ldiag[Nx] = 0
        ldiag[-1] = 0
        ldiag = ldiag[1:]  # discard the 0th entry

        # A lonely upper diagonal to incorporate G_B, the BC for A depending on b
        ulonely = np.zeros(Nx)
        ulonely[0] = -G_B(t)

        # Two entries to be added to account the conservation of mass (BC of b)
        # For the moment I could only do it by creating two diagonals
        first = np.zeros(Nx)
        first[0] = -1
        second = np.zeros(Nx + 1)
        second[1] = 1

        # This creates a sparse matrix out of specified diagonals
        return diags(
            [diag, udiag, ulonely, ldiag, first, second],
            [0, 1, Nx, -1, -Nx, -Nx + 1],
            format="csr",
        )

    ## Solve the unknowns iteratively with time
    # First create an empty vector for each quantity a and b and set the first values to be the IC
    sol_A = np.zeros(Nx * Nt)
    sol_A[0:Nx] = IC_A

    sol_B = np.zeros(Nx * Nt)
    sol_B[0:Nx] = IC_B

    # Loop over all time-step
    for t in range(1, Nt):
        matrix = oneTimeStepMatrix(t)
        print(matrix)
        rhs = np.concatenate(
            ([0], sol_A[Nx*(t - 1) + 1 : Nx * (t - 1) + Nx - 1], [1], [0], sol_B[Nx*(t - 1) + 1 : Nx*(t - 1) + Nx - 1], [0]))
        x = spsolve(matrix, rhs)
        sol_A[Nx*t : Nx*t + Nx] = x[0:Nx]
        sol_B[Nx*t : Nx*t + Nx] = x[Nx:]
    plotAnimate(x0, xn, Nx, Nt, dt, sol_A, sol_B, "volt.gif","Concentration A","Concentration B")
    return sol_A, sol_B


def plotAnimate(x0, xn, Nx, Nt, dt, x_A, x_B, filename, x_A_Str,x_B_Str):
    # Plotting
    fig, ax = plt.subplots(figsize=(5, 3))
    minTot = min(x_B.min(),x_A.min())
    maxTot = max(x_B.max(),x_A.max())
    ax.set(xlim=(x0, xn), ylim=(minTot, maxTot))
    x = np.linspace(x0, xn, Nx)

    line = ax.plot(x, x_A[0:Nx], color="k", lw=2,label = x_A_Str)[0]  # x_A(0), ...x_A(Nx-1)
    line2 = ax.plot(x, x_B[0:Nx], color="r", lw=2,label = x_B_Str)[0]  # x_A(0), ...x_A(Nx-1)

    plt.legend(loc="lower right")
    plt.title("Analytical versus Numerical Solution")
    plt.xlabel("thtnh")
    plt.ylabel("concentration")

    def animate(t):
        first = t * Nx
        last = first + Nx - 1
        # print("Frame", x_A[first : last + 1].mean()) #who wroe that?
        line.set_ydata(x_A[first : last + 1])
        line2.set_ydata(x_B[first : last + 1])

    anim = FuncAnimation(fig, animate, interval=dt * 400, frames=Nt)
    anim.save(str(RESULTS_FOLDER / filename))
    plt.show()
    return x_A, x_B


if __name__ == "__main__":
    # chronoamperometry()
    [x_A,x_B]=voltametry()
    print(x_A[0:10])
    print(x_B[0:10])

