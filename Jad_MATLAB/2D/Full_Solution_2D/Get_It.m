function [It] = Get_It (Z,dx,dy,ay,x1)
%This function computes the current by integrating the normal derivative
%of concentration a(U variable) along y=0 only for x<x1
l = length(Z(1,round(ay/dx)+1:round(x1/dx)+1));
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
