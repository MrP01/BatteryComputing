function [ZU,ZV] = create_Z(U_plot,V_plot,Nx,Ny,c,BCx2U,BCy2U,BCx2V,BCy2V,ICU,dx,Uy0_U,Uy0_V,Ux0_U,Ux0_V,x,Upy0_U,Upy0_V,Vy0_U,Vy0_V,Vx0_U,Vx0_V,Vpy0_U,Vpy0_V,ICV)
%This function creates both ZU and ZV which are the values of U and V to be 
% drawn on the grid 

%Defining constants to be used:
    x1 = round(x/dx);
    l1 = length(Uy0_U);
    l2 = length(Uy0_V);
    l3 = length(Ux0_U);
    l4 = length(Ux0_V);
    l5 = length(Upy0_U);
    l6 = length(Upy0_V);
    ZU = zeros(Ny+2,Nx+2);
    ZV = zeros(Ny+2,Nx+2);

%Specifying the Dirichlet BCs 
    ZU(end,:)=ones(1,Nx+2)*BCy2U;
    ZV(end,:)=ones(1,Nx+2)*BCy2V;
    ZU(:,end) = ones(Ny+2,1)*BCx2U;
    ZV(:,end) = ones(Ny+2,1)*BCx2V;

    counter = 0;
    if c~=0
    %Filling Inner points
    for i = 2:Ny+1
    ZU(i,2:Nx+1)=U_plot(1+(Nx*Ny)*(c-1)+Nx*counter:Nx+(Nx*Ny)*(c-1)+Nx*counter,1);
    ZV(i,2:Nx+1)=V_plot(1+(Nx*Ny)*(c-1)+Nx*counter:Nx+(Nx*Ny)*(c-1)+Nx*counter,1);
    counter = counter+1;
    end
    %The following function fills values for BC on x=0 and y=0
    [ZU]=fillZ(Nx,Ny,Ux0_U,Ux0_V,Upy0_U,Upy0_V,Uy0_U,Uy0_V,x1,l1,l2,l3,l4,l5,l6,ZU,ZV,1);
    l1 = length(Vy0_U);
    l2 = length(Vy0_V);
    l3 = length(Vx0_U);
    l4 = length(Vx0_V);
    l5 = length(Vpy0_U);
    l6 = length(Vpy0_V);
    %The following function fills values for BC on x=0 and y=0
    [ZV]=fillZ(Nx,Ny,Vx0_U,Vx0_V,Vpy0_U,Vpy0_V,Vy0_U,Vy0_V,x1,l1,l2,l3,l4,l5,l6,ZU,ZV,0);
    else
    %Only for c=0, the values are inititialized to be the initial
    %conditions
        ZU = ones(size(ZU))*ICU;
        ZV = ones(size(ZV))*ICV;
    end
end