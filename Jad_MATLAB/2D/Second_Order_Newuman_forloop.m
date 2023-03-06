clc
clearvars
tic
ax = 0;
bx = 10;
ay = 0;
by = 10;
T0 = 0;
T = 3;
dx = 0.05;
dy = 0.05;
dt = 0.1;
NBCx1 = 0;
NBCy1 = 0;
BCx2 = 1;
BCy2 = 1;
BCy1 = 0;
IC = 1;
[dx,dy,dt,Nx,Ny,Nt] = create_mesh(ax,bx,ay,by,T0,T,dx,dy,dt);
U = ones(1,Nx*Ny)*IC;
U_plot = zeros(Nx*Ny*Nt,1);
counter = 0;
[A] = Create_matrix_f(Nx,Ny,dt/(dx)^2,dt/(dy)^2,dx);
for t = T0:dt:T-dt
[bs] = Create_RHS_f(Nx,Ny,dt/(dx)^2,dt/(dy)^2,NBCx1,NBCy1,BCx2,BCy1,BCy2,dx,dy,dt,U);
U = (A\bs)';
U_plot(1+Nx*Ny*counter:Nx*Ny*(counter+1),1) = U;
counter = counter+1;
clc
fprintf('Working on it: %f percent',t/(T-dt)*100)
end
fprintf('\n')
[X,Y]=meshgrid(ax:dx:bx,ay:dy:by);
[Z] = create_Z(U_plot,Nx,Ny,0,NBCx1,BCx2,NBCy1,BCy1,BCy2,IC,dx,dy);
toc
figure
surf(X,Y,Z)
xlim([ax,bx])
ylim([ay,by])
zlim([min(U_plot)*(1+sign(min(U_plot))*0.2),IC*1.3])
counter = 1;
I = zeros(1,Nt);
for t = T0:dt:T-dt
    pause(0.1)
    [Z] = create_Z(U_plot,Nx,Ny,counter,NBCx1,BCx2,NBCy1,BCy1,BCy2,IC,dx,dy);
    I(1,counter)=Get_It(Z,dx,dy);
    counter = counter+1;
    surf(X,Y,Z)
    hold on
    xlim([ax,bx])
    ylim([ay,by])
    zlim([min(U_plot)*(1+sign(min(U_plot))*0.2),IC*1.3])
    hold off
end

figure
plot(T0:dt:T-dt,I)
legend('Current')
title('Current in function of time')
xlabel('Time')
ylabel('Current')
%% Functions
function [It] = Get_It (Z,dx,dy)
l = length(Z(1,:));
    I = zeros(1,l);
    for i = 1:l
        I(1,i) = (-3/2*Z(1,i)+2*Z(2,i)-0.5*Z(3,i))/dy;
    end
sum = 0;
    for i = 2:l
        sum = sum+(I(1,i)+I(1,i-1))/2*dx;
    end
It = sum;
end


function [dx,dy,dt,Nx,Ny,Nt] = create_mesh(ax,bx,ay,by,T0,T,dx,dy,dt)
    Nx = round(abs((ax-bx))/dx)-1;
    Ny = round(abs((ay-by))/dy)-1;
    Nt = round(abs((T-T0))/dt);
    dx = abs((ax-bx))/(Nx+1);
    dy = abs((ay-by))/(Ny+1);
    dt = abs((T-T0))/(Nt);
end

function [A] = Create_matrix_f (Nx,Ny,mx,my,dx)
x1 = round(1/dx);
r1 = 1:Nx*Ny;
c1 = r1;
nnz1 = ones(1,length(r1))*(1+mx+my);
rd = nonzeros((1:Nx*Ny).*(rem((1:Nx*Ny)-1,Nx)~=0))';
r2 = [rd,rd];
c2 = [rd,rd-1];
nnz2 = [ones(1,length(rd))*mx,ones(1,length(rd))*-mx];
r3 = nonzeros((1:Nx*Ny).*(rem((1:Nx*Ny),Nx)~=0))';
c3 = r3+1;
nnz3 = ones(1,length(r3))*-mx;
rd = Nx+1:Nx*Ny;
r4 = [rd,rd];
c4 = [rd,rd-Nx];
nnz4 = [ones(1,length(rd))*my,ones(1,length(rd))*-my];
r5 = 1:(Nx-1)*Ny;
c5 = r5+Nx;
nnz5 = ones(1,length(r5))*-my;
r6 = 1:x1;
c6 = r6;
nnz6 = ones(1,length(r6))*my;
r = [r1,r2,r3,r4,r5,r6];
c = [c1,c2,c3,c4,c5,c6];
nnz = [nnz1,nnz2,nnz3,nnz4,nnz5,nnz6];
A = sparse(r,c,nnz);
end

function [b] = Create_RHS_f(Nx,Ny,mx,my,NBCx1,NBCy1,BCx2,BCy1,BCy2,dx,dy,dt,U)
x1 = round(1/dx);
r1 = 1:Nx*Ny;
r2 = nonzeros((1:Nx*Ny).*(rem((1:Nx*Ny)-1,Nx)==0))';
r3 = nonzeros((1:Nx*Ny).*(rem((1:Nx*Ny),Nx)==0))';
r4 = x1+1:Nx;
r5 = (Nx-1)*Ny+1:Nx*Ny;
r6 = 1:x1;

r = [r1,r2,r3,r4,r5,r6];
c = ones(1,length(r));
nnz =[U,ones(1,length(r2))*-dt/dx*NBCx1,ones(1,length(r3))*BCx2*mx,ones(1,length(r4))*-dt/dy*NBCy1,ones(1,length(r5))*my*BCy2,ones(1,length(r6))*my*BCy1];
b = sparse(r,c,nnz);
end

function [Z] = create_Z(U,Nx,Ny,c,NBCx1,BCx2,NBCy1,BCy1,BCy2,IC,dx,dy)
    x1 = round(1/dx);
    Z = zeros(Ny+2,Nx+2);
    Z(end,:)=ones(1,Nx+2)*BCy2;
    Z(end,1)= (NBCx1*dx+0.5*BCy2-2*BCy2)*-2/3;
    Z(end,1) = -(NBCx1*dx-BCy2);
    Z(:,end) = ones(Ny+2,1)*BCx2;
    Z(1,1:x1) = ones(1,x1)*BCy1;
    counter = 0;
    if c~=0
    for i = 2:Ny+1
    Z(i,2:Nx+1)=U(1+(Nx*Ny)*(c-1)+Nx*counter:Nx+(Nx*Ny)*(c-1)+Nx*counter,1);
    counter = counter+1;
    end
    for i = 2:Ny+1
    Z(i,1)= (NBCx1*dx+0.5*Z(i,2)-2*Z(i,3))*-2/3;
    end
    for i = x1+1:Nx+1
    Z(1,i) = (NBCy1*dy+0.5*Z(2,i)-2*Z(3,i))*-2/3;
    end
    else
        Z = ones(size(Z))*IC;
    end
end
