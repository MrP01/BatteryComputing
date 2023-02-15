function[rA1,cA1,nnzA1,rA2,cA2,nnzA2] = Sparsity_for_BCx0yo(Uy0_U,Uy0_V,Vy0_U,Vy0_V,Upy0_U,Upy0_V,Vpy0_U,Vpy0_V,Ux0_U,Ux0_V,Vx0_U,Vx0_V,Nx,Ny,x,ay,by,ax,bx,dx,dy,D,mx,my)

%This function returns the entries of the matrix related to the BC
%specified on x=0 an y=0:

[rA1a,cA1a,nnzAU1a,nnzAV1a,~,~,~,~] = create_sparsity_pattern(Uy0_U,Uy0_V,Vy0_U,Vy0_V,Nx,Ny,"y=0",x+dx,by-dx,dx);

if x ~=ay %for x<x1
[rA1b,cA1b,nnzAU1b,nnzAV1b,~,~,~,~] = create_sparsity_pattern(Upy0_U,Upy0_V,Vpy0_U,Vpy0_V,Nx,Ny,"y=0",ay+dx,x,dx);
else
rA1b = []; 
cA1b = [];
nnzAU1b = [];
nnzAV1b = [];
end

rA1 = [rA1a,rA1b];
cA1 = [cA1a,cA1b];
nnzAU1 = [nnzAU1a,nnzAU1b];
nnzAV1 = [nnzAV1a,nnzAV1b];
[rA2,cA2,nnzAU2,nnzAV2,~,~,~,~] = create_sparsity_pattern(Ux0_U,Ux0_V,Vx0_U,Vx0_V,Nx,Ny,"x=0",ax+dy,bx-dy,dy);
nnzA1 = [-my*nnzAU1,-D*my*nnzAV1];
nnzA2 = [-mx*nnzAU2,-D*mx*nnzAV2];
end