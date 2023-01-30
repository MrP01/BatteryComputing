# Python code for solving chronoamperometre problem

from scipy.sparse import dia_matrix, dok_matrix
from scipy.sparse.linalg import spsolve
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

## Set up for solving the heat equation for a

# Boundaries for space and time
x0 = 0 ; xn = 10 # i replaced a by x_0 and b by xn
T0 = 0 ; T = 1

# Space mesh
dx = 0.1 ; dt = 0.01

# Values at boundary in space and time
BC1 = 0; BC2 = 1
IC = 1

# Set up matrices for solving the scheme
[dx,dt,Nx,Nt] = create_mesh(a,b,T0,T,dx,dt);
[A] = Create_matrix_Dirichlet(Nx-1,Nt,dt/(dx)^2);
[bs] = Create_RHS_Dirichlet(Nx-1,Nt,dt/(dx)^2,BC1,BC2,IC);
U1 = A\bs;
U1 = full(U1);
counter = 0;
