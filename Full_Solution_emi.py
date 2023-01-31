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
#dx = 0.1 ; dt = 0.01

# Number of meshpoints and meshsizes
Nx=100; Nt=100
dx=(xn-x0)/Nx ; dt = (T-T0)/Nt

# Courant Friedrichs-Lewy
CFL = dt/dx^2;

# Values at boundary in space and time
BC1 = 0; BC2 = 1
IC = 1

# Set up matrices for solving the scheme
[A] = schemeDirichlet(Nx-1,Nt,CFL);
[bs] = rhsDirichlet(Nx-1,Nt, CFL,BC1,BC2,IC);
U1 = A\bs;
U1 = full(U1);
counter = 0;

# Function creating the scheme matrix for a
def schemeDirichlet(Nx,Nt,m)
    diag=np.zeros((Nx+1)*(Nt+1))+1
    lonelydiag=np.zeros((Nx+1)*(Nt+1)-(Nx+1))+1
    upperdiag=np.zeros((Nx+1)*(Nt+1)-1)
    lowerdiag = np.zeros((Nx + 1) * (Nt + 1) - 1)
    for t in range(1,Nt+1)  #so for times t=1,2,...Nt
        lonelydiag(t*Nx+1)=0
        lonelydiag((t+1)*Nx-1)=0
        for k in range(2, Nx) #so for k=2...Nx-1
            diag(t*Nx+k) = 1+2*CFL
            upperdiag(t*Nx+k) = -CFL
            lowerdiag(t*Nx+k-1)=-CFL

    for

    %% Everything except the first row
    rl = nonzeros([Nx+1:Nx*Nt].*(rem([Nx+1:Nx*Nt]-1,Nx)~=0))';
    cl = rl-1;
    rr = nonzeros([Nx+1:Nx*Nt].*(rem([Nx+1:Nx*Nt],Nx)~=0))';
    cr = rr+1;

r1 = [Nx+1:Nx*Nt,Nx+1:Nx*Nt,rl,rr];
c1 = [1:Nx*(Nt-1),(Nx+1:Nx*Nt),cl,cr];

l1 = length(Nx+1:Nx*Nt);
l2 = length(Nx+1:Nx*Nt);
l3 = length(rl);
l4 = length(rr);
nnz1 = [-ones(1,l1),ones(1,l2)*(1+2*m),-ones(1,l3)*m,-ones(1,l4)*m];
%% First row
r0 = [1:Nx,1:Nx-1,2:Nx];
c0 = [1:Nx,(1:Nx-1)+1,(2:Nx)-1];
nnz0 = [ones(1,Nx)*(1+2*m),-ones(1,Nx-1)*m,-ones(1,Nx-1)*m];

%% Create the matrix
A = sparse([r1,r0],[c1,c0],[nnz1,nnz0]);

return A

