function [rUA,cUA,nnzUA,rUb,nnzUb] = sparsity (vecs,veco,Nx,Ny,s1,s2,x0,xf,dx)

if s1 == "y=0"
       x0 = round(x0/dx);
       xf = round(xf/dx);
       N = xf-x0+1;
    if s2 == "U"
    rd = x0:xf;
    elseif s2 == "V"
    rd = (x0:xf)+(Nx*Ny);
    end
elseif s1 == "x=0"
            N = Ny;
    rd = nonzeros((1:Nx*Ny).*(rem((1:Nx*Ny)-1,Nx)==0))';
    if s2 == "V"
    rd = rd+Nx*Ny;
    end
else
        error('S value can only be defined as \"x=0\" or \"y=0\"')
end
l = length(vecs)-1;
ru1 = zeros(1,l*N);
cu1 = zeros(1,l*N);
nnzu1 = zeros(1,l*N);
for i = 1:l
    ru1(1,1+N*(i-1):N+N*(i-1)) = rd;
    cu1(1,1+N*(i-1):N+N*(i-1)) = rd+(i-1);
    nnzu1(1,1+N*(i-1):N+N*(i-1)) = vecs(1,i+1);
end
if s2 == "U"
    Inc = Nx*Ny;
elseif s2 == "V"
    Inc = -Nx*Ny;
end
nnzu2 = zeros(1,(length(veco)-1)*N);
for i = 1:length(veco)-1
    ru2(1,1+N*(i-1):N+N*(i-1)) = rd;
    cu2(1,1+N*(i-1):N+N*(i-1)) = rd+(i-1)+Inc;
    nnzu2(1,1+N*(i-1):N+N*(i-1)) = veco(1,i+1);
end
if isempty(nnzu2)
    ru2 = [];
    cu2 = [];
end

rUA = [ru1,ru2];
cUA = [cu1,cu2];
nnzUA = [nnzu1,nnzu2];
rUb = rd;
nnzUb = ones(1,length(rUb))*(vecs(1,1)+veco(1,1));
end