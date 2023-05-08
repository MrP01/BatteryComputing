function [ZZ]=fillZ(Nx,Ny,Zx0_U,Zx0_V,Zpy0_U,Zpy0_V,Zy0_U,Zy0_V,x1,l1,l2,l3,l4,l5,l6,ZU,ZV,c)
if c ==1
    ZZ = ZU;
else
    ZZ = ZV;
end
    for i = 2:Ny+1
        sum = Zx0_U(1,1);
        for n = 2:l3
            sum = sum+Zx0_U(1,n)*ZU(i,n);
        end
        sum = sum+Zx0_V(1,1);
        for n = 2:l4
            sum = sum+Zx0_V(1,n)*ZV(i,n);
        end
    ZZ(i,1)= sum;
    end
    for i = 1:Nx+1
        if i<=x1
    sum = Zpy0_U(1,1);
        for n = 2:l5
            sum = sum+Zpy0_U(1,n)*ZU(n,i);
        end
        sum = sum+Zpy0_V(1,1);
        for n = 2:l6
            sum = sum+Zpy0_V(1,n)*ZV(n,i);
        end
    ZZ(1,i) = sum;
        else
        sum = Zy0_U(1,1);
        for n = 2:l1
            sum = sum+Zy0_U(1,n)*ZU(n,i);
        end
        sum = sum+Zy0_V(1,1);
        for n = 2:l2
            sum = sum+Zy0_V(1,n)*ZV(n,i);
        end
    ZZ(1,i) = sum;
        end
    end
        ZZ(1,1) = (ZZ(2,1)+ZZ(1,2))/2;
end
