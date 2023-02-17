function [dx,dy,dt,Nx,Ny,Nt] = create_mesh(ax,bx,ay,by,T0,T,dx,dy,dt)
%This function creates the mesh for certain bounds and increments and gives
%back the number of internal points and exact increments.
    Nx = round(abs((ax-bx))/dx)-1;
    Ny = round(abs((ay-by))/dy)-1;
    Nt = round(abs((T-T0))/dt);
    dx = abs((ax-bx))/(Nx+1);
    dy = abs((ay-by))/(Ny+1);
    dt = abs((T-T0))/(Nt);
end
