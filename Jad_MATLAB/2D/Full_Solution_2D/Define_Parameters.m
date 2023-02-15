function [ax,bx,ay,by,T0,T,dx,dy,dt,D,Uy0_U,Uy0_V,Ux0_U,Ux0_V,Vy0_U,Vy0_V,Vx0_U,Vx0_V,BCx2U,BCy2U,BCx2V,BCy2V,ICU,ICV,x1] = Define_Parameters()
%In this function, every parameter is defined with brief comment explaining
%what they are
ax = 0; %Lowest x value for U and V
bx = 10; %Highest x value for U and V 
ay = ax;
by = bx;
T0 = 0; %Starting time
T = 5; %Ending time
dx = 0.5; %Space increment
dy = dx; 
dt = 0.1; %Time increment
D = 1; %Check explanation (diffusion coeff)
x1 = 1.5;%Specifying where BC at y=0 change

%NOTE: U represent a and V represent b

Uy0_U = [0 1]; %U0 as a function of [cte, U1, U2,...] at y = 0 (for x>=x1)
Uy0_V = 0; %U0 as a function of [cte, V1, V2,...] at y = 0 (for x>=x1)
Ux0_U = [0 1]; %U0 as a function of [cte, U1, U2,...] at x = 0
Ux0_V = 0; %U0 as a function of [cte, V1, V2,...] at y = 0

%The following are the same for V
Vy0_U = 0;
Vy0_V = [0 1];
Vx0_U = 0;
Vx0_V = [0 1];

%Specifying BC at higher x and ys:
BCx2U = 1;
BCy2U = 1;
BCx2V = 0;
BCy2V = 0;

%Specifying initial conditions:
ICU = 1;
ICV = 0;
end