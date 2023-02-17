function [rA,cA,nnzAU,nnzAV,rbU,rbV,nnzbU,nnzbV] = create_sparsity_pattern(UvecU,VvecU,UvecV,VvecV,Nx,Ny,s,x0,xf,dx)
    [rA1,cA1,nnzAU,rbU,nnzbU] = sparsity (UvecU,VvecU,Nx,Ny,s,"U",x0,xf,dx);
    [rA2,cA2,nnzAV,rbV,nnzbV] = sparsity (VvecV,UvecV,Nx,Ny,s,"V",x0,xf,dx);
    rA = [rA1,rA2];
    cA = [cA1,cA2];
end