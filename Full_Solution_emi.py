# Python code for solving chronoamperometre problem

import math
import pathlib
import numpy as np

import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from scipy import special
from scipy.sparse import diags
from scipy.sparse.linalg import spsolve

# Latex rendering
plt.rcParams['text.usetex'] = True
plt.rcParams["font.family"] = "mathpazo"

# Python code for External viewer on Mac
# plt.switch_backend("MacOSX")
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
        first = t * Nx  # first boundary entry of a timestep
        last = first + Nx - 1  # last boundary entry of a timestep

        diag[first] = 1  # first BC entry of a timestep
        diag[last] = 1  # last BC entry of a timestep

        ldiag[last] = 0  # last boundary entry
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
def rhsNeumann(Nx, Nt, IC, BC2, a, D):  # Neumann for the first boundary and Dirichlet with the second
    r = np.zeros(Nx * Nt)
    # Initial conditions for the first Nx rows of r.
    r[0:Nx] = IC

    # Fill the boundary conditions in
    for t in range(1, Nt):
        r[t * Nx] = -D * (3 / 2 * a[t * Nx] - 2 * a[t * Nx + 1] + 1 / 2 * a[t * Nx + 2])  # using the Taylor
        r[t * Nx + Nx - 1] = BC2
    return r


## CHRONOAMPEROMETRY is the main code for running experiment with Dirichlet BC for A and Neumann for B
def chronoamperometry(x0, xn, t0, tn, Nx, dx, Nt, dt, D):
    # Courant Friedrichs-Lewy
    CFL = dt / dx**2

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
    real_A = np.zeros(Nx * Nt)
    real_A[0:Nx] = IC_A
    for j in range(1, Nt):
        x = np.linspace(0, xn, Nx)
        h = math.sqrt(j * dt)
        real_A[j * Nx : (j + 1) * Nx] = special.erf(x / (2 * h))

    # Plot only one instance
"""     plotChrono(x0, xn, Nx, Nt, dt, x_A, x_B, "chrono_AvsB_", "Numerical Solution", "Real Solution",0, False)
    plotChrono(x0, xn, Nx, Nt, dt, x_A, x_B, "chrono_AvsB_", "Numerical Solution", "Real Solution",1,False)
    plotChrono(x0, xn, Nx, Nt, dt, x_A, x_B, "chrono_AvsB_", "$A$", "$B$", Nt-1,True)
  """
    # Compare real to numerical solution
    #plotAnimate(x0, xn, Nx, Nt, dt, x_A, real_A, "chrono_AvsReal.gif", "Numerical Solution", "Real Solution")

    # Plotting the concentration of A versus that of B
    #plotAnimate(x0, xn, Nx, Nt, dt, x_A, x_B, "chrono_AvsB.gif", "Concentration A", "Concentration B")


# VOLTAMETRY is the main code for running our second experiment, namely using the Butler-Volmer equation
def voltametry(x0, xn, t0, tn, Nx, dx, Nt, dt, D, alpha, DC, delta, ome):
    # Courant Friedrichs-Lewy and diffusion constant
    CFL = dt / dx**2
    D = 1

    # Values at boundary in space and time
    IC_A = 1
    IC_B = 0

    # Initial values
    E_start = -10
    E_0 = 0
    t_rev = 20
    kappa = 35

    def Pot(t, DC):  # This function returns our potential
        E = E_start + t if t < t_rev else E_start - t + 2 * t_rev
        if not DC:
            E = E + delta * math.sin(ome * t)
        return E

    def G_A(t, DC):  # This function represents the first term of our Butler-Volmer function
        return kappa * dx * math.exp((1 - alpha) * (Pot(t, DC) - E_0))

    def G_B(t, DC):  # This function represents the positive part of the second term of our Butler-Volmer function
        return kappa * dx * math.exp((-alpha) * (Pot(t, DC) - E_0))

    # oneTimeStepMatrix creates a matrix that incorporates the Butler-Volmer scheme but only for one step.
    # This is similar to schemeD, only that schemeD creates a matrix involving all unknowns in time and space
    def oneTimeStepMatrix(j):
        ## Create a matrix for one timestep. The unknowns are (A0...A_N-1,B0,...B_N-1)
        # The diagonal is 1+2*D*mu except at the boundaries
        diag = np.ones(2 * Nx) * (1 + 2 * D * CFL)
        diag[0] = G_A(dt * j, DC) + 1
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
        ulonely[0] = -G_B(j * dt, DC)

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

    # Then create arrays for the current I and the potential E
    I_AC = np.zeros(Nt)  # I[0]=0
    I = np.zeros(Nt)
    E = np.zeros(Nt)
    E[0] = Pot(0, True)  # Storing the potential values for DC

    # Loop over all time-step
    for j in range(1, Nt):
        matrix = oneTimeStepMatrix(j)
        rhs = np.concatenate(
            (
                [0],
                sol_A[Nx * (j - 1) + 1 : Nx * (j - 1) + Nx - 1],
                [1],
                [0],
                sol_B[Nx * (j - 1) + 1 : Nx * (j - 1) + Nx - 1],
                [0],
            )
        )
        x = spsolve(matrix, rhs)
        sol_A[Nx * j : Nx * j + Nx] = x[0:Nx]
        sol_B[Nx * j : Nx * j + Nx] = x[Nx:]
        # Compute the current
        I[j] = (x[1] - x[0]) / dx
        E[j] = Pot(j * dt, True)

    return sol_A, sol_B, I, E


def plotAnimate(x0, xn, Nx, Nt, dt, x_A, x_B, filename, x_A_Str, x_B_Str):
    # Plotting
    fig, ax = plt.subplots(figsize=(5, 5))
    minTot = min(x_B.min(), x_A.min())
    maxTot = max(x_B.max(), x_A.max())
    ax.set(xlim=(x0, xn), ylim=(minTot, maxTot))
    x = np.linspace(x0, xn, Nx)

    line = ax.plot(x, x_A[0:Nx], color="k", lw=2, label=x_A_Str)[0]  # x_A(0), ...x_A(Nx-1)
    line2 = ax.plot(x, x_B[0:Nx], color="r", lw=2, label=x_B_Str)[0]  # x_A(0), ...x_A(Nx-1)

    plt.legend(loc="lower right")
    plt.xlabel("x")
    plt.ylabel("concentration")

    def animate(j):
        first = j * Nx
        last = first + Nx - 1
        line.set_ydata(x_A[first : last + 1])
        line2.set_ydata(x_B[first : last + 1])
        nb = j * dt
        plt.title("Analytical versus Numerical Solution at t=" + str("%.2f" % nb))

    anim = FuncAnimation(fig, animate, interval=dt * 400, frames=Nt)
    anim.save(str(RESULTS_FOLDER / filename))
    # plt.show()
    return x_A, x_B

def plotChrono(x0, xn, Nx, Nt, dt, x_A, x_B, filetype, x_A_Str, x_B_Str, timeFrame_j,leg):
    # Plotting
    fig, ax = plt.subplots(figsize=(5, 5))
    minTot = min(x_B.min(), x_A.min())
    maxTot = max(x_B.max(), x_A.max())
    ax.set(xlim=(x0, xn), ylim=(minTot, maxTot))
    x = np.linspace(x0, xn, Nx)

    line = ax.plot(x, x_A[0:Nx], color="k", lw=2, label=x_A_Str)[0]  # x_A(0), ...x_A(Nx-1)
    line2 = ax.plot(x, x_B[0:Nx], color="r", lw=2, label=x_B_Str)[0]  # x_A(0), ...x_A(Nx-1)

    if leg:
        plt.legend(loc="lower right", fontsize=20)

    plt.xlabel("$x$", fontsize = 20)
    plt.ylabel("Concentration", fontsize = 20)

    j=timeFrame_j

    first = j * Nx
    last = first + Nx - 1
    line.set_ydata(x_A[first : last + 1])
    line2.set_ydata(x_B[first : last + 1])
    

    # Generate (time, the time float) (timeStr, the string from time) (intTime , 100 times time) (intTimeStr, string for intTime)
    time = j * dt
    timeStr = str("%.2f" % time)

    intTime = int(time*100)
    intTimeStr = str(intTime)
    plt.title("$t=$" + timeStr, fontsize = 25)

    # Save some figures
    filename = filetype + intTimeStr + ".png"
    plt.savefig(str(RESULTS_FOLDER / filename))

    # plt.show()
    return x_A, x_B


# Used for plotting
def myPlotting(x0, xn, Nx, Nt, dt, D, AC, delta, omega):
    # for j in range(0,11):
    # alpha = j *0.1
    alpha = 0.2
    [x_A, x_B, I, E] = voltametry(x0, xn, t0, tn, Nx, dx, Nt, dt, D, alpha, AC, delta, omega)
    plt.figure()
    plt.plot(E, I, label=alpha)
    plt.xlabel("E")
    plt.ylabel("I")
    plt.title("alpha=" + str(alpha) + ", T_rev=" + str(20) + ", T=" + str(tn))
    plt.savefig(str(RESULTS_FOLDER / ("volt_alpha" + str(alpha) + ".png")))


if __name__ == "__main__":
    ## Set up for solving the heat equation for a
    # Boundaries for space and time
    x0 = 0
    xn = 20
    t0 = 0
    tn = 10

    # Number of meshpoints and meshsizes
    Nx = 500
    dx = (xn - x0) / (Nx - 1)
    Nt = 600
    dt = (tn - t0) / (Nt - 1)

    # Other parameters
    alpha = 0.2
    D = 1
    DC = True
    delta = 0.1
    omega = 2 * math.pi

    #chronoamperometry(x0, xn, t0, tn, Nx, dx, Nt, dt, D)
    [sol_A, sol_B, I, E] = voltametry(x0, xn, t0, tn, Nx, dx, Nt, dt, D, alpha, DC, delta, omega)

    # Plotting the concentration of A versus that of B
     if DC:
         myStr = "volt_DC_AvsB.gif"
     else:
         myStr = "volt_AC_AvsB.gif"

     plotAnimate(x0, xn, Nx, Nt, dt, sol_A, sol_B, myStr, "Concentration A", "Concentration B")

    # # Plotting for different alphas
     myPlotting(x0, xn, Nx, Nt, dt, D, DC, delta, omega)
