
function [b] = Create_RHS_f(Nx,Ny,mx,my,BCx2U,BCy2U,BCx2V,BCy2V,dx,dy,U,V,D,Uy0_U,Uy0_V,Ux0_U,Ux0_V,Vy0_U,Vy0_V,Vx0_U,Vx0_V,ax,bx,ay,by,x,Upy0_U,Upy0_V,Vpy0_U,Vpy0_V)

%This function does the same thing as creating A but for the right hand
%side vector b!

[~,~,~,~,rbU1a,rbV1a,nnzbU1a,nnzbV1a] = create_sparsity_pattern(Uy0_U,Uy0_V,Vy0_U,Vy0_V,Nx,Ny,"y=0",x+dx,by-dx,dx);
if x~=ay
[~,~,~,~,rbU1b,rbV1b,nnzbU1b,nnzbV1b] = create_sparsity_pattern(Upy0_U,Upy0_V,Vpy0_U,Vpy0_V,Nx,Ny,"y=0",ay+dx,x,dx);
else
rbU1a = [];
rbV1a = [];
nnzbU1a = [];
nnzbV1a = [];
rbU1b = [];
rbV1b = [];
nnzbU1b = [];
nnzbV1b = [];
end
[~,~,~,~,rb2a,rb2b,nnzbU2,nnzbV2] = create_sparsity_pattern(Ux0_U,Ux0_V,Vx0_U,Vx0_V,Nx,Ny,"x=0",ax+dy,bx-dy,dy);
rb2 = [rb2a,rb2b];
r1U = nonzeros((1:Nx*Ny).*(rem((1:Nx*Ny),Nx)==0))';
r2U = (Nx-1)*Ny+1:Nx*Ny;
r3U = 1:Nx*Ny;

r1V = r1U+Nx*Ny;
r2V = r2U+Nx*Ny; 
r3V = r3U+Nx*Ny;

r = [r1U,r2U,r3U,r1V,r2V,r3V,rbU1a,rbV1a,rbU1b,rbV1b,rb2];
c = ones(1,length(r));
nnz =[ones(1,length(r1U))*BCx2U*mx,ones(1,length(r2U))*my*BCy2U,U,ones(1,length(r1V))*BCx2V*mx*D,ones(1,length(r2V))*my*BCy2V*D,V,my*nnzbU1a,my*D*nnzbV1a,my*nnzbU1b,my*D*nnzbV1b,mx*nnzbU2,mx*D*nnzbV2]; 
b = sparse(r,c,nnz);
end