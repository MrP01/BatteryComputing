function [I] = Plot_Results(Draw,U_plot,V_plot,Nx,Ny,BCx2U,BCy2U,BCx2V,BCy2V,ICU,dx,Uy0_U,Uy0_V,Ux0_U,Ux0_V,x1,Upy0_U,Upy0_V,Vy0_U,Vy0_V,Vx0_U,Vx0_V,Vpy0_U,Vpy0_V,ICV,ax,bx,ay,by,dy,dt,Nt,T0,T)
[X,Y]=meshgrid(ax:dx:bx,ay:dy:by);
[ZU,ZV] = create_Z(U_plot,V_plot,Nx,Ny,0,BCx2U,BCy2U,BCx2V,BCy2V,ICU,dx,Uy0_U,Uy0_V,Ux0_U,Ux0_V,x1,Upy0_U,Upy0_V,Vy0_U,Vy0_V,Vx0_U,Vx0_V,Vpy0_U,Vpy0_V,ICV);

toc
if Draw ~=0
figure(1) %Drawing U
clf
surf(X,Y,ZU)
xlim([ax,bx])
ylim([ay,by])
zlim([0,max(max(ZU))+0.1])
figure(2) %Drawing V
clf
surf(X,Y,ZV)
xlim([ax,bx])
ylim([ay,by])
zlim([0,max(max(ZV)+0.1)])
figure(3) %Drawing U+V
surf(X,Y,ZU+ZV)
xlim([ax,bx])
ylim([ay,by])
zlim([0,max(max(ZU+ZV))+0.1])
end
counter = 1;
I = zeros(1,Nt);

for t = T0:dt:T-dt
    
    
    [ZU,ZV] = create_Z(U_plot,V_plot,Nx,Ny,counter,BCx2U,BCy2U,BCx2V,BCy2V,ICU,dx,Uy0_U,Uy0_V,Ux0_U,Ux0_V,x1,Upy0_U,Upy0_V,Vy0_U,Vy0_V,Vx0_U,Vx0_V,Vpy0_U,Vpy0_V,ICV);
    I(1,counter)=Get_It(ZU,dx,dy,ay,x1);
    counter = counter+1;
    if Draw~=0
    pause(0.1)
    figure(1)
    surf(X,Y,ZU)
    hold on
    xlim([ax,bx])
    ylim([ay,by])
    zlim([0,max(max(ZU))+0.1])
    hold off
    figure(2)
    surf(X,Y,ZV)
    hold on
    xlim([ax,bx])
    ylim([ay,by])
    zlim([0,max(max(ZV))+0.1])
    hold off
    figure(3)
    surf(X,Y,ZU+ZV)
    hold on
    xlim([ax,bx])
    ylim([ay,by])
    zlim([0,max(max(ZU+ZV))+0.1])
    hold off
    end
end
if Draw ~=0
figure
plot(T0:dt:T-dt,I)
legend('Current')
title('Current in function of time')
xlabel('Time')
ylabel('Current')
end
end