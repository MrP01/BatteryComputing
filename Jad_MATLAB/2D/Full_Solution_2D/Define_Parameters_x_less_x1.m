function [Upy0_U,Upy0_V,Vpy0_U,Vpy0_V]= Define_Parameters_x_less_x1 (c,D,t,T,dy)
%This function specifies parameters in vector form (just like in define
%parameters previously) for x<x1
%Data = [E0 k0 alpha w]

    if c~=1
        Estart = -10;
        E0 = 0;
        trev = T/2;
        k0 = 35;
        alpha = 0.5;%Changed from 0.5
        DE = 0.1;
        w = 2*pi;
    if t<=trev
        E = @(t) Estart+t+DE*sin(w*t);
    else
        E = @(t) Estart+trev-(t-trev)+DE*sin(w*t);
    end
        f1 = @(t) k0*exp((1-alpha)*(E(t)-E0));
        f2 = @(t) k0*exp(-alpha*(E(t)-E0));
        C = D*(1+dy*f1(t))+dy*f2(t);
        Upy0_U = [0,1-dy*f1(t)*D];
        Upy0_V = [0,D-D*(1+dy*f1(t))*D];
        Vpy0_U = [0,dy*f1(t)/C];
        Vpy0_V = [0,D*(1+dy*f1(t))/C];
    else
        Upy0_U = 0;
        Upy0_V = 0;
        Vpy0_U = [0,1/D];
        Vpy0_V = [0 1];
    end
end
