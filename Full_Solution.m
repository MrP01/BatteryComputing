clc;
clearvars;
a = 0;
b = 10;
T0 = 0;
T = 1;
dx = 0.1;
dt = 0.01;
BC2 = 1;
IC = 1;
BC1 = 0;
[dx,dt,Nx,Nt] = create_mesh(a,b,T0,T,dx,dt);
[A] = Create_matrix_Dirichlet(Nx-1,Nt,dt/(dx)^2);
[bs] = Create_RHS_Dirichlet(Nx-1,Nt,dt/(dx)^2,BC1,BC2,IC);
U1 = A\bs;
U1 = full(U1);
counter = 0;
% figure 
% plot(a:dx:b,ones(1,Nx+1)*IC)
% xlim([a,b])
% ylim([0,IC*1.3])
D = zeros(1,length(T0:dt:T-2*dt));
counter = 0;
for t = T0+dt:dt:T
%     pause(0.1)
%     plot(a:dx:b,[BC1,U1(1+(Nx-1)*counter:(Nx-1)+(Nx-1)*counter)',BC2])
%     hold on
%     if t~=0
%     plot(a:dx:b,erf((a:dx:b)/(2*sqrt(t))))
%     end
%     hold off  
    D(1,counter+1)= ((-3/2)*BC1+2*U1(1+(Nx-1)*counter)-0.5*U1(2+(Nx-1)*counter))/dx;
%     xlim([a,b])
%     ylim([0,IC*1.3])
%     if t~=0
%     legend('Approx','Solution')
%     end
    counter =counter+1;
end
Dp = [((-3/2)*IC+2*IC-0.5*IC)/dx,D];
% figure
% plot(T0:dt:T-dt,D)
%hold on
%plot(T0:dt:T-dt,1-D)
%%

BC2n = 0;
IC = 0;
Cu = 1;
NBC1 = -D/Cu;
[A] = Create_matrix_Newmann(Nx-1,Nt,dt/(dx)^2);
[bs] = Create_RHS_Newmann(Nx-1,Nt,dt/(dx)^2,NBC1,BC2n,IC,dt,dx);
U = A\bs;
U = full(U);
counter = 0;
figure 
% plot(a:dx:b,ones(1,Nx+1)*IC)
% xlim([a,b])
% ylim([0,1.3])
%D = zeros(1,length(T0:dt:T-dt));
for t = T0+dt:dt:T
    pause(0.1)
    plot(a:dx:b,[-(NBC1(1,counter+1)*dx-U(1+(Nx-1)*counter)),U(1+(Nx-1)*counter:(Nx-1)+(Nx-1)*counter)',BC2n],'LineWidth',2)
    hold on
    plot(a:dx:b,[BC1,U1(1+(Nx-1)*counter:(Nx-1)+(Nx-1)*counter)',BC2],'LineWidth',2)
    if Cu == 1
    plot(a:dx:b,erf((a:dx:b)/(2*sqrt(t))),'LineWidth',2)
    plot(a:dx:b,1-erf((a:dx:b)/(2*sqrt(t))),'LineWidth',2)
    end
    hold off
    xlim([a,b])
    ylim([0,1.3])
    if Cu == 1
    legend('b','a','exact_a','exact_b')
    else
    legend('b','a')
    end
    counter =counter+1;
end
% D =[(-3/2)*BC1+2*IC-0.5*IC,D];
figure
plot(T0:dt:T,Dp)


function [dx,dt,Nx,Nt] = create_mesh(a,b,T0,T,dx,dt)
    Nx = round(abs((a-b))/dx);
    Nt = round(abs((T-T0))/dt);
    dx = abs((a-b))/(Nx);
    dt = abs((T-T0))/(Nt);
end

function [A] = Create_matrix_Dirichlet (Nx,Nt,m)
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

function [A] = Create_matrix_Newmann (Nx,Nt,m)
%% Everything except the first row
rl1 = nonzeros([Nx+1:Nx*Nt].*(rem([Nx+1:Nx*Nt]-1,Nx)~=0).*(rem([Nx+1:Nx*Nt],Nx)~=0))';

r1 = [Nx+1:Nx*Nt,rl1,rl1,rl1];
c1 = [1:Nx*(Nt-1),rl1,rl1-1,rl1+1];
l1 = length(Nx+1:Nx*Nt);
l3 = length(rl1);
nnz1 = [-ones(1,l1),ones(1,l3)*(1+2*m),-ones(1,l3)*m,-ones(1,l3)*m];

rl2 = nonzeros([Nx+1:Nx*Nt].*(rem([Nx+1:Nx*Nt]-1,Nx)==0))';
r2 = [rl2,rl2];
c2 = [rl2,rl2+1];
l3 = length(rl2);
nnz2 = [ones(1,l3)*(1+m),-ones(1,l3)*m];

rl3 = nonzeros([Nx+1:Nx*Nt].*(rem([Nx+1:Nx*Nt],Nx)==0))';
r3 = [rl3,rl3];
c3 = [rl3,rl3-1];
l3 = length(rl3);
nnz3 = [ones(1,l3)*(1+2*m),-ones(1,l3)*m];
%% First row
r0 = [1,1,2:Nx-1,2:Nx-1,2:Nx-1,Nx,Nx];
c0 = [1,2,2:Nx-1,(2:Nx-1)+1,(2:Nx-1)-1,Nx-1,Nx];
nnz0 = [1+m,-m,ones(1,Nx-2)*(1+2*m),-ones(1,Nx-2)*m,-ones(1,Nx-2)*m,-m,1+2*m];

%% Create the matrix
A = sparse([r1,r2,r3,r0],[c1,c2,c3,c0],[nnz1,nnz2,nnz3,nnz0]);     
end


function [b] = Create_RHS_Dirichlet(Nx,Nt,m,BC1,BC2,IC)
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

function [b] = Create_RHS_Newmann(Nx,Nt,m,NBC1,BC2,IC,dt,dx)

rl = nonzeros([1:Nx*Nt].*(rem([1:Nx*Nt]-1,Nx)==0))';
cl = ones(1,length(rl));
nnzl = -dt/dx*NBC1;
rr = nonzeros([1:Nx*Nt].*(rem([1:Nx*Nt],Nx)==0))';
cr = ones(1,length(rr));
nnzr =  ones(1,length(rr))*m*BC2;
ro = 1:Nx;
co = ones(1,Nx);
nnzo = ones(1,Nx)*IC;

b = sparse([rl,rr,ro],[cl,cr,co],[nnzl,nnzr,nnzo]);
end
