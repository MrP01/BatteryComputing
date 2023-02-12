clc
clearvars
tic
ax = 0;
bx = 10;
ay = 0;
by = 10;
T0 = 0;
T = 1;
dx = 0.05;
dy = 0.05;
dt = 0.1;
NBCx1 = 0.1;
NBCy1 = 0.1;
BCx2 = 0;
BCy2 = 0;
IC = 1;
[dx,dy,dt,Nx,Ny,Nt] = create_mesh(ax,bx,ay,by,T0,T,dx,dy,dt);
[A] = Create_matrix(Nx,Ny,Nt,dt/(dx)^2,dt/(dy)^2);
[bs] = Create_RHS(Nx,Ny,Nt,dt/(dx)^2,dt/(dy)^2,NBCx1,NBCy1,BCx2,BCy2,IC,dx,dy,dt);
U = A\bs;
U = full(U);
[X,Y]=meshgrid(ax:dx:bx,ay:dy:by);
[Z] = create_Z(U,Nx,Ny,0,NBCx1,BCx2,NBCy1,BCy2,IC,dx,dy);
toc
figure 
surf(X,Y,Z)
xlim([ax,bx])
ylim([ay,by])
zlim([min(U)*0.8,IC*1.3])
counter = 1;
for t = T0:dt:T-dt
    pause(0.1)
    [Z] = create_Z(U,Nx,Ny,counter,NBCx1,BCx2,NBCy1,BCy2,IC,dx,dy);
    counter = counter+1;
    surf(X,Y,Z)
    hold on
    xlim([ax,bx])
    ylim([ay,by])
    zlim([0,IC*1.3])
    hold off
end


function [dx,dy,dt,Nx,Ny,Nt] = create_mesh(ax,bx,ay,by,T0,T,dx,dy,dt)
    Nx = round(abs((ax-bx))/dx)-1;
    Ny = round(abs((ay-by))/dy)-1;
    Nt = round(abs((T-T0))/dt);
    dx = abs((ax-bx))/(Nx+1);
    dy = abs((ay-by))/(Ny+1);
    dt = abs((T-T0))/(Nt);
end

function [A] = Create_matrix (Nx,Ny,Nt,mx,my)
r1 = 1:Nx*Ny*Nt;
c1 = r1;
nnz1 = ones(1,length(r1))*(1+mx+my);
rd = nonzeros([1:Nx*Ny*Nt].*(rem([1:Nx*Ny*Nt]-1,Nx)~=0))';
r2 = [rd,rd];
c2 = [rd,rd-1];
nnz2 = [ones(1,length(rd))*mx,ones(1,length(rd))*-mx];
r3 = nonzeros([1:Nx*Ny*Nt].*(rem([1:Nx*Ny*Nt],Nx)~=0))';
c3 = r3+1;
nnz3 = ones(1,length(r3))*-mx;
r0 = (rem([1:Nx*Ny*Nt]-1,Nx*Ny)~=0);
for i = 2:Nx
    r0=r0.*(rem([1:Nx*Ny*Nt]-i,Nx*Ny)~=0);
end
rd = nonzeros([1:Nx*Ny*Nt].*r0)';
r4 = [rd,rd];
c4 = [rd,rd-Nx];
nnz4 = [ones(1,length(rd))*my,ones(1,length(rd))*-my];
r0 = (rem([1:Nx*Ny*Nt]+-((Nx-1)*Ny+1),Nx*Ny)~=0);
for i = (Nx-1)*Ny+2:Nx*Ny
    r0=r0.*(rem([1:Nx*Ny*Nt]-i,Nx*Ny)~=0);
end
r5 = nonzeros([1:Nx*Ny*Nt].*r0)';
c5 = r5+Nx;
nnz5 = ones(1,length(r5))*-my;
r6 = Nx*Ny+1:Nx*Ny*Nt;
c6 = 1:length(r6);
nnz6 = ones(1,length(r6))*-1;
r = [r1,r2,r3,r4,r5,r6];
c = [c1,c2,c3,c4,c5,c6];
nnz = [nnz1,nnz2,nnz3,nnz4,nnz5,nnz6];
A = sparse(r,c,nnz);     
end

function [b] = Create_RHS(Nx,Ny,Nt,mx,my,NBCx1,NBCy1,BCx2,BCy2,IC,dx,dy,dt)

r2 = nonzeros([1:Nx*Ny*Nt].*(rem([1:Nx*Ny*Nt]-1,Nx)==0))';
r3 = nonzeros([1:Nx*Ny*Nt].*(rem([1:Nx*Ny*Nt],Nx)==0))';
r0 = (rem([1:Nx*Ny*Nt]-1,Nx*Ny)==0);
for i = 2:Nx
    r0=r0+(rem([1:Nx*Ny*Nt]-i,Nx*Ny)==0);
end
r4 = nonzeros([1:Nx*Ny*Nt].*(r0~=0))';
r0 = (rem([1:Nx*Ny*Nt]-((Nx-1)*Ny+1),Nx*Ny)==0);
for i = (Nx-1)*Ny+2:Nx*Ny
    r0=r0+(rem([1:Nx*Ny*Nt]-i,Nx*Ny)==0);
end
r5 = nonzeros([1:Nx*Ny*Nt].*(r0~=0))';
r6 = 1:Nx*Ny;

r = [r2,r3,r4,r5,r6];
c = ones(1,length(r));
nnz =[ones(1,length(r2))*-dt/dx*NBCx1,ones(1,length(r3))*BCx2*mx,ones(1,length(r4))*-dt/dy*NBCy1,ones(1,length(r5))*my*BCy2,ones(1,length(r6))*IC]; 
b = sparse(r,c,nnz);
end

function [Z] = create_Z(U,Nx,Ny,c,NBCx1,BCx2,NBCy1,BCy2,IC,dx,dy)
    Z = zeros(Ny+2,Nx+2);
    Z(end,:)=ones(1,Nx+2)*BCy2;
    Z(end,1)= (NBCx1*dx+0.5*BCy2-2*BCy2)*-2/3;
    Z(end,1) = -(NBCx1*dx-BCy2);
    Z(:,end) = ones(Ny+2,1)*BCx2;
    counter = 0;
    if c~=0
    for i = 2:Ny+1
    Z(i,2:Nx+1)=U(1+(Nx*Ny)*(c-1)+Nx*counter:Nx+(Nx*Ny)*(c-1)+Nx*counter,1);
    counter = counter+1;
    end
    for i = 2:Ny+1
    Z(i,1)= (NBCx1*dx+0.5*Z(i,2)-2*Z(i,3))*-2/3;
    Z(i,1) = -(NBCx1*dx-Z(i,2));
    end
    for i = 2:Nx+1
    Z(1,i) = (NBCy1*dy+0.5*Z(2,i)-2*Z(3,i))*-2/3;
    Z(1,i) = -(NBCy1*dy-Z(2,i));
    end
       Z(1,1) = (Z(1,2)+Z(2,1))/2;
    else
        Z = ones(size(Z))*IC;
    end
end

