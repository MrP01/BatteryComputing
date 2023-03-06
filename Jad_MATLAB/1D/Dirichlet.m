clc;

a = 0;
b = 10;
T0 = 0;
T = 3;
dx = 0.1;
dt = 0.01;
BC2 = 1;
IC = 1;
BC1 = 0;
[dx,dt,Nx,Nt] = create_mesh(a,b,T0,T,dx,dt);
[A] = Create_matrix(Nx-1,Nt,dt/(dx)^2);
[bs] = Create_RHS(Nx-1,Nt,dt/(dx)^2,BC1,BC2,IC);
U = A\bs;
U = full(U);
counter = 0;
figure
plot(a:dx:b,ones(1,Nx+1)*IC)
xlim([a,b])
ylim([0,IC*1.3])
D = zeros(1,length(T0:dt:T-dt));
for t = T0:dt:T-dt
    pause(0.1)
    plot(a:dx:b,[BC1,U(1+(Nx-1)*counter:(Nx-1)+(Nx-1)*counter)',BC2])
    hold on
    if t~=0
    plot(a:dx:b,erf((a:dx:b)/(2*sqrt(t))))
    end
    hold off
    D(1,counter+1)= (-3/2)*BC1+2*U(1+(Nx-1)*counter)-0.5*U(2+(Nx-1)*counter);
    xlim([a,b])
    ylim([0,IC*1.3])
    if t~=0
    legend('Approx','Solution')
    end
    counter =counter+1;
end
D =[(-3/2)*BC1+2*IC-0.5*IC,D];
figure
plot(T0:dt:T,D)
hold on
plot(T0:dt:T,1-D)
x = a+dx:dx:b-dx;
for t = T0+dt:dt:T-dt
    x = [x,[a+dx:dx:b-dx]];
end



function [dx,dt,Nx,Nt] = create_mesh(a,b,T0,T,dx,dt)
    Nx = round(abs((a-b))/dx);
    Nt = round(abs((T-T0))/dt);
    dx = abs((a-b))/(Nx);
    dt = abs((T-T0))/(Nt);
end

function [A] = Create_matrix (Nx,Nt,m)
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
end

function [b] = Create_RHS(Nx,Nt,m,BC1,BC2,IC)
rl = nonzeros([Nx+1:Nx*Nt].*(rem([Nx+1:Nx*Nt]-1,Nx)==0))';
cl = ones(1,length(rl));
nnzl = ones(1,length(rl))*m*BC1;
rr = nonzeros([Nx+1:Nx*Nt].*(rem([Nx+1:Nx*Nt],Nx)==0))';
cr = ones(1,length(rr));
nnzr =  ones(1,length(rr))*m*BC2;
ro = 1:Nx;
co = ones(1,Nx);
nnzo = [m*BC1+IC,ones(1,Nx-2)*IC,m*BC2+IC];

b = sparse([rl,rr,ro],[cl,cr,co],[nnzl,nnzr,nnzo]);
end
