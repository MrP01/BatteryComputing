function [Upy0_U,Upy0_V,Vpy0_U,Vpy0_V]= Define_Parameters_x_less_x1 (c,D)
%This function specifies parameters in vector form (just like in define
%parameters previously) for x<x1
    if c~=1
    if t<=T/2
        E = @(t) -10+t;
    else
        E = @(t) -10+T/2-(t-T/2);
    end
        f1 = @(t) 35*exp((1-0.5)*(E(t)));
        f2 = @(t) 35*exp(-0.5*(E(t)));
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