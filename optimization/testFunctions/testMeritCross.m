function [ e ] = testMeritCross( x,v,u,obj,ss,varargin)


opt = struct('pert',1e-5,'debug',false);
opt = merge_options(opt, varargin{:});


withAlgs = ss.nv>0;

pertV = opt.pert;

e = [];

xDims = cellfun(@(zi)numel(zi),x);

f = rand;
dE = {cellfun(@(it)rand(size(it))-0.5,x,'UniformOutput',false)};
rho = 1;



[ m,Jm ] = l1merit(f,dE,rho,'gradients',true);

mf = @(ff) l1merit(ff,dE,rho);
dmdf = calcPertGrad(mf,f,pertV);

lag = @(dd) l1merit(f,{mat2cell(dd,xDims,1)},rho);
dmdd = calcPertGrad(lag,cell2mat(dE{1}),pertV);


e = [ e norm(Jm.Jf-dmdf)];
e = [ e norm(cell2mat(Jm.JdE)-dmdd)];

fSeed = rand;
dSeed = cellfun(@(it)rand(size(it)),x,'UniformOutput',false);
ls = rand;


[ m,JmR ] = l1merit(f,dE,rho,'gradients',true,'fRightSeeds',fSeed,'dERightSeeds',{dSeed},'leftSeed',ls);

JmRF =  ls*(Jm.Jf*fSeed+ cell2mat(Jm.JdE)*cell2mat(dSeed));
e = [e norm(JmRF-JmR.J)];


rho = 1;


    
merit = @(f,dE,varargin) l1merit(f,dE,rho,varargin{:});





[xs,vs,xd,vd,a,A,av,Av] = condensing(x,u,v,ss,'computeCorrection',true);


du = cellfun(@(z)rand(size(z))/10,u,'UniformOutput',false);
dxN = cellmtimes(A,du,'lowerTriangular',true,'ci',ss.ci);
if withAlgs
    dvN = cellmtimes(Av,du,'lowerTriangular',true,'ci',ss.ci);
else
    dvN = [];
end
dx = cellfun(@(z,dz)z+dz,a,dxN,'UniformOutput',false);
if withAlgs
    dv = cellfun(@(z,dz)z+dz,av,dvN,'UniformOutput',false);
else
    dv =[];
end


simFunc = @(xk,uk,varargin) simulateSystem(xk,uk,ss,varargin{:});
phi = @(l,varargin) lineFunctionWrapper(l,x,v,u,dx,dv,du,...
        simFunc,obj,merit,'gradients',true,'plotFunc',[],'plot',false,...
        varargin{:});


[f0,g0]=  lineFunctionWrapper(0,x,v,u,dx,dv,du,...
        simFunc,obj,merit,'gradients',true,'plotFunc',[],'plot',false,...
        'xd0',xd,'vd0',vd,'xs0',xs,'vs0',vs);


e = [e,norm((phi(opt.pert)-f0)/opt.pert-g0)/norm(g0)];

ll = rand;

[f,g]=  phi(ll);
e = [e,norm((phi(opt.pert+ll)-f)/opt.pert-g)/norm(g)];





end



function [dfdx] = calcPertGrad(f,xb,pert)

nx = numel(xb);
fb = f(xb);

dfdx = zeros(numel(fb),nx);

parfor j = 1:nx
    xp = xb;
    xp(j) = xp(j) + pert;
    
    fp = f(xp);
    
    dfdx(:,j) = (fp-fb)/pert;
    
end

end

