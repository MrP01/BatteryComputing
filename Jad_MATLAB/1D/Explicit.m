clc;
clearvars;

a = 0;
b = 1;
T0 = 0;
T = 10;
dx = 10^-2;
dt = 0.5*10^-4;
BC2 = 0;
IC = 1;
BC1 = 0;
[U,dx,dt] = create_mesh(a,b,T0,T,dx,dt,IC);
counter = 0;

figure
for t = T0:dt:T-dt
    U = march_in_time(U,dx,dt,BC1,BC2);
    plot (a:dx:b,U)
    xlim([a,b])
    ylim([0,IC*1.1])
    pause(0.1)
end



function [U,dx,dt] = create_mesh(a,b,T0,T,dx,dt,IC)
    Nx = round(abs((a-b))/dx)+1;
    Nt = round(abs((T-T0))/dt)+1;
    dx = abs((a-b))/(Nx-1);
    dt = abs((T-T0))/(Nt-1);
    U = ones(1,Nx)*IC;
end

function [U] = march_in_time(U,dx,dt,BC1,BC2)
        U = U +dt*(diffusion(U,dx));
        U(1,1) = BC1;
        U(1,end) = BC2;
end

function [d] = diffusion(U,dx)
d = zeros(1,length(U));
    for x = 2:length(U)-1
        d(1,x) = (U(1,x-1)-2*U(1,x)+U(1,x+1))/dx^2;
    end
end