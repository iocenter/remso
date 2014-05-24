function [ maxerror ] = testLiftOptReduction(x,u,v,ss )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

nx = numel(x{1});
nu = numel(u{1});
nv = ss.nv;
withAlgs  = nv >0;



[xsF,vF,Jac] = simulateSystem(x,u,ss,'gradients',true);

d = cellfun(@minus,xsF,x,'UniformOutput',false);

dJac = Jac.xJx;
for k = 1:numel(x)
    dJac{k,k} = -eye(numel(x{k}));
end

% lift-opt related matrix brute force calculated
invdhmidx = inv(cell2matFill(dJac,[nx,nx]));
a = -invdhmidx*cell2mat(d);
A  = -invdhmidx*cell2matFill(Jac.xJu,[nx,nu]);
%b = fg + cell2mat(fgx)*a;

[xs,vs,xd,vd,axZ,asu,av,Av] = condensing(x,u,v,ss);




e1 = norm(cell2mat(axZ)-a);
e2 = norm(full(cell2matFill(asu,[nx,nu]))-A);
e3 = norm(cell2mat(d)-cell2mat(xd));
if nv >0
    e4 =   norm(cell2mat(av) - (cell2mat(vd) + cell2matFill(Jac.vJx,[nv,nx])*a));
    e5 =   norm(full(cell2matFill(Av,[nv,nu]) - (cell2matFill(Jac.vJu,[nv,nu])+cell2matFill(Jac.vJx,[nv,nx])*A)));
else
    e4 = 0;
    e5 = 0;
end
if withAlgs
    e6 =   norm(cell2mat(vs)-cell2mat(v)-cell2mat(vd));
else
    e6 =0;
end

e7 = norm(cell2mat(xsF)-cell2mat(xs));
if withAlgs
    e8 = norm(cell2mat(vs) - cell2mat(vF));
else
    e8 = 0;
end
    
    
maxerror = max([e1,e2,e3,e4,e5,e6,e7,e8]);

end

