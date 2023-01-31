# Python code for solving chronoamperometre problem

from scipy.sparse import dia_matrix, dok_matrix
from scipy.sparse.linalg import spsolve
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

## Set up for solving the heat equation for a

# Boundaries for space and time
x0 = 0 ; xn = 1 # i replaced a by x_0 and b by xn #should be xn=10
T0 = 0 ; T = 1

# Space mesh
#dx = 0.1 ; dt = 0.01

# Number of meshpoints and meshsizes
Nx=10; Nt=10
dx=(xn-x0)/Nx ; dt = (T-T0)/Nt

# Courant Friedrichs-Lewy
CFL = dt/dx^2;

# Values at boundary in space and time
BC1 = 0; BC2 = 1
IC = 1

# Set up matrices for solving the scheme
[A] = schemeDirichlet(Nx-1,Nt,CFL);
#[bs] = rhsDirichlet(Nx-1,Nt, CFL,BC1,BC2,IC);
#U1 = A\bs;
#U1 = full(U1);
#counter = 0;

# Function creating the scheme matrix for a
def schemeDirichlet(Nx,Nt,m)

    #diagonal vector with (Nx+1)*(Nt+1) entries
    diag=np.zeros((Nx+1)*(Nt+1))+1

    #lonely diagonal vector accounts for the explicit coefficient
    lonelydiag=np.zeros((Nx+1)*(Nt+1)-(Nx+1))+1

    #off diagonal vectors
    udiag=np.zeros((Nx+1)*(Nt+1)-1)
    ldiag = np.zeros((Nx + 1) * (Nt + 1) - 1)

    # We aim to zero the entries of the lonely vector corresponding to boundary
    # We aim to have tridiagonal matrix (lowerdiag, diag, upperdiag) with (-CFl, 1+2CFL,CFL)
    for t in range(1,Nt+1)  #so for times t=1,2,...Nt
        lonelydiag(t*Nx+1)=0
        lonelydiag((t+1)*Nx-1)=0
        for k in range(2, Nx) #so for k=2...Nx-1
            diag(t*Nx+k) = 1+2*CFL
            udiag(t*Nx+k) = -CFL
            ldiag(t*Nx+k-1)=-CFL

    A= dia_matrix(([diag, udiag, ldiag, lonelydiag], [0, 1, -1, Nx+1]),shape=(Nx+1, Nt+1))
return A

print(A.toarray())