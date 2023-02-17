function [A] = Create_matrix_f (Nx,Ny,mx,my,dx,dy,D,Uy0_U,Uy0_V,Ux0_U,Ux0_V,Vy0_U,Vy0_V,Vx0_U,Vx0_V,x,ax,bx,ay,by,Upy0_U,Upy0_V,Vpy0_U,Vpy0_V)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%The following function returns a matrix with entries for all a values and
%all b values at a specific time step to be solved (implicit scheme)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[rA1,cA1,nnzA1,rA2,cA2,nnzA2] = Sparsity_for_BCx0yo(Uy0_U,Uy0_V,Vy0_U,Vy0_V,Upy0_U,Upy0_V,Vpy0_U,Vpy0_V,Ux0_U,Ux0_V,Vx0_U,Vx0_V,Nx,Ny,x,ay,by,ax,bx,dx,dy,D,mx,my);

%Entries for U(i,j)
rU1 = 1:Nx*Ny;
cU1 = rU1;
nnzU1 = ones(1,Nx*Ny)*(1+2*mx+2*my);

%Entries for U(i,j-1)
rU2 = Nx+1:Nx*Ny;
cU2 = rU2-Nx;
nnzU2 = ones(1,length(rU2))*(-my);

%Entries for U(i-1,j)
rU3 = nonzeros((1:Nx*Ny).*(rem((1:Nx*Ny)-1,Nx)~=0))';
cU3 = rU3-1;
nnzU3 = ones(1,length(rU3))*(-mx);

%Entries for U(i,j+1)
rU4 = 1:(Nx-1)*Ny;
cU4 = rU4+Nx;
nnzU4 = ones(1,length(rU4))*(-my);

%Entries for column next to boundary x=L
rU5 = nonzeros((1:Nx*Ny).*(rem((1:Nx*Ny),Nx)~=0))';
cU5 = rU5+1;
nnzU5 = ones(1,length(rU5))*(-mx);

%Entries for row just before x=L
rV1 = Nx*Ny+1:2*Nx*Ny;
cV1 = rV1;
nnzV1 = ones(1,Nx*Ny)*(1+D*(2*mx+2*my));

%Entries for V(i,j-1)
rV2 = rU2+Nx*Ny;
cV2 = rV2-Nx;
nnzV2 = ones(1,length(rV2))*(-D*my);

%Entries for V(i-1,j)
rV3 = rU3+Nx*Ny;
cV3 = rV3-1;
nnzV3 = ones(1,length(rV3))*(-D*mx);

%Entries for V(i,j+1)
rV4 = rU4+Nx*Ny;
cV4 = rV4+Nx;
nnzV4 = ones(1,length(rV4))*(-D*my);

%Entries for V(i+1,j)
rV5 = rU5+Nx*Ny;
cV5 = rV5+1;
nnzV5 = ones(1,length(rV5))*(-D*mx);

r = [rA1,rA2,rU1,rU2,rU3,rU4,rU5,rV1,rV2,rV3,rV4,rV5];
c = [cA1,cA2,cU1,cU2,cU3,cU4,cU5,cV1,cV2,cV3,cV4,cV5];
nnz = [nnzA1,nnzA2,nnzU1,nnzU2,nnzU3,nnzU4,nnzU5,nnzV1,nnzV2,nnzV3,nnzV4,nnzV5];

A = sparse(r,c,nnz);     
end