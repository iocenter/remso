function [ maxError ] = testNonlinearGradient(x,u,v,ss,obj,varargin)


opt = struct('pert',1e-5,'debug',false);
opt = merge_options(opt, varargin{:});

maxError = -inf;
pertV = opt.pert;

nx = numel(x{1});
nu = numel(u{1});
nv = numel(v{1});

withAlgs = (ss.nv >0);


%x = cellfun(@(xi)xi+norm(xi)*rand(size(xi))/20,x,'UniformOutput',false);
%v = cellfun(@(xi)xi+norm(xi)*rand(size(xi))/20,v,'UniformOutput',false);
%u = cellfun(@(xi)xi+norm(xi)*rand(size(xi))/20,u,'UniformOutput',false);



[xs,vs,xd,vd,a,A,av,Av] = condensing(x,u,v,ss);


[f,gradU] = targetGrad(x,u,v,obj,A,Av,ss.ci);


%{
%% this is being tested elsewhere

[xsF,vF,JacFull] = simulateSystem(x,u,ss,'gradients',true);

xJu = cell2matFill(JacFull.xJu,[nx,nu]);
xJx = cell2matFill(JacFull.xJx,[nx,nx]);
vJu = cell2matFill(JacFull.vJu,[nv,nu]);
vJx = cell2matFill(JacFull.vJx,[nv,nx]);

xd = cell2mat(cellfun(@minus,xsF,x,'UniformOutput',false));

invX = inv(xJx-eye(size(xJx) ));
a = -invX*xd;
A  = -invX*xJu;


dxdu = (-xJx * invX  + eye(size(xJx)))*xJu;
dvdu = vJu + vJx*dxdu;


[f,targetPartials] = obj(x,u,v,'gradients',true);

dfdu = cell2mat(targetPartials.Jx)*dxdu + cell2mat(targetPartials.Jv)*dvdu + cell2mat(targetPartials.Ju);


norm(cell2mat(gradU)-dfdu)
norm(A-dxdu)

%}




[fz,Bz,cz,simz,zz,zv,usz] = simulateSystemZ(u,xd,vd,ss,obj,'gradients',true);



e = [];
e = [e f-fz norm(cell2mat(gradU)-cell2mat(Bz))];




objPartials.Jx = cellfun(@(xi)rand([size(xi,2),size(xi,1)]),x','UniformOutput',false);
objPartials.Jv = cellfun(@(xi)rand([size(xi,2),size(xi,1)]),v','UniformOutput',false);
objPartials.Ju = cellfun(@(xi)rand([size(xi,2),size(xi,1)]),u','UniformOutput',false);
%[f,objPartials] = obj(x,u,v,'gradients',true);

if withAlgs
    gZ = vectorTimesZ(objPartials.Jx,objPartials.Ju,objPartials.Jv,A,Av,ss.ci );
else
    gZ = vectorTimesZ(objPartials.Jx,objPartials.Ju,[],A,[],ss.ci );
end


[xsR,vsR,~,convergedR,simVarsR,uslicedR] = simulateSystem(x,u,ss,'gradients',false,'guessX',xs,'guessV',vs);
xdR = cellfun(@(xsi,xi)xsi-xi,xsR,x,'UniformOutput',false);
vdR = cellfun(@(xsi,xi)xsi-xi,vsR,v,'UniformOutput',false);



[~,gradUY,converged,~,xs,vs,zusliced ] = simulateSystemZ(u,xdR,vdR,ss,[],'gradients',true,'guessV',xsR,'guessX',vsR,'simVars',simVarsR,'JacTar',objPartials);


w = cellfun(@minus,gradUY,gZ,'UniformOutput',false);


e = [e norm(cell2mat(w))];




uSeed = cellfun(@(it)rand(size(it)),u,'UniformOutput',false);



[~,~,~,~,~,AxS,~,AvS] = condensing(x,u,v,ss,'uRightSeeds',uSeed);


e = [e, norm(cell2matFill(A,[nx,nu])*cell2mat(uSeed)-cell2mat(AxS)),norm(cell2matFill(Av,[nv,nu])*cell2mat(uSeed)-cell2mat(AvS)) ];




[f,JacTarFull] = obj(xs,u,vs,'gradients',true);


fx = @(xit) obj(toStructuredCells(xit,nx),u,vs);
dfdx = calcPertGrad(fx,cell2mat(xs),pertV);


fu = @(uit) obj(xs,toStructuredCells(uit,nu),vs);
dfdu = calcPertGrad(fu,cell2mat(u),pertV);

fv = @(vit) obj(xs,u,toStructuredCells(vit,nv));
if withAlgs
    dfdv = calcPertGrad(fv,cell2mat(vs),pertV);
end

e = [e norm(cell2mat(JacTarFull.Jx)-dfdx)];
e = [e norm(cell2mat(JacTarFull.Ju)-dfdu)];

if withAlgs
    e = [e norm(cell2mat(JacTarFull.Jv)-dfdv)];
end

xSeed = cellfun(@(it)rand(size(it)),xs,'UniformOutput',false);
uSeed = cellfun(@(it)rand(size(it)),u,'UniformOutput',false);
vSeed = cellfun(@(it)rand(size(it)),vs,'UniformOutput',false);
leftSeed = rand(size(f))';


[f,JacTarRight] = obj(xs,u,vs,'gradients',true,'xRightSeeds',xSeed,'uRightSeeds',uSeed,'vRightSeeds',vSeed);
fullJ = cell2mat(JacTarFull.Jx)*cell2mat(xSeed)+cell2mat(JacTarFull.Ju)*cell2mat(uSeed)+cell2mat(JacTarFull.Jv)*cell2mat(vSeed);
e = [e norm(fullJ-JacTarRight.J)];

[f,JacTarLeft] = obj(xs,u,vs,'gradients',true,'leftSeed',leftSeed);
e = [e norm(leftSeed*cell2mat(JacTarFull.Jx)-cell2mat(JacTarLeft.Jx))];
e = [e norm(leftSeed*cell2mat(JacTarFull.Ju)-cell2mat(JacTarLeft.Ju))];
e = [e norm(leftSeed*cell2mat(JacTarFull.Jv)-cell2mat(JacTarLeft.Jv))];

[f,gradU] = targetGrad(xs,u,vs,obj,A,Av,ss.ci);   %%% x or xs for this evaluation??  --> test

if withAlgs
    gradUF = cell2mat(JacTarFull.Ju) + cell2mat(JacTarFull.Jx) * cell2matFill(A,[nx,nu]) + cell2mat(JacTarFull.Jv) * cell2matFill(Av,[nv,nu]);
else
    gradUF = cell2mat(JacTarFull.Ju) + cell2mat(JacTarFull.Jx) * cell2matFill(A,[nx,nu]) ;
end
e = [e norm(cell2mat(gradU)-gradUF)];

%
% %
% %
% %
% fobj = @(uit) msObjective(x,toStructuredCells(uit,nu),ss,obj);
% dfdu = calcPertGrad(fobj,cell2mat(u),pertV);
% %
% dfdu-cell2mat(gradU)
%
% %finitErr = 0;
%
% finitErr = norm(dfdu-cell2mat(gradU));


f = rand;
dE = {cellfun(@(it)rand(size(it))-0.5,x,'UniformOutput',false)};
bE = {cellfun(@(it)rand(size(it)),x,'UniformOutput',false)};
lb = {cellfun(@(it)rand(size(it))-0.5,x,'UniformOutput',false)};
ub = {cellfun(@(it)rand(size(it))+0.5,x,'UniformOutput',false)};
rho = 1;
tau = {cellfun(@(it)rand(size(it)),x,'UniformOutput',false)};



[ m,Jm ] = l1merit(f,dE,bE,ub,lb,rho,tau,'gradients',true);

mf = @(ff) l1merit(ff,dE,bE,ub,lb,rho,tau);
dmdf = calcPertGrad(mf,f,pertV);

md = @(dd) l1merit(f,{toStructuredCells(dd,nx)},bE,ub,lb,rho,tau);
dmdd = calcPertGrad(md,cell2mat(dE{1}),pertV);

mb = @(bb) l1merit(f,dE,{toStructuredCells(bb,nx)},ub,lb,rho,tau);
dmdb = calcPertGrad(mb,cell2mat(bE{1}),pertV);

e = [ e norm(Jm.Jf-dmdf)];
e = [ e norm(cell2mat(Jm.JdE)-dmdd)];
e = [ e norm(cell2mat(Jm.JbE)-dmdb)];

fSeed = rand;
dSeed = cellfun(@(it)rand(size(it)),x,'UniformOutput',false);
bSeed = cellfun(@(it)rand(size(it)),x,'UniformOutput',false);
ls = rand;


[ m,JmR ] = l1merit(f,dE,bE,ub,lb,rho,tau,'gradients',true,'fRightSeeds',fSeed,'dERightSeeds',{dSeed},'bERightSeeds',{bSeed},'leftSeed',ls);

JmRF =  ls*(Jm.Jf*fSeed+ cell2mat(Jm.JdE)*cell2mat(dSeed)+ cell2mat(Jm.JbE)*cell2mat(bSeed));
e = [e norm(JmRF-JmR.J)];


rho = 1;

if withAlgs
    lb = [{cellfun(@(it)rand(size(it))-0.5,x,'UniformOutput',false)};{cellfun(@(it)rand(size(it))-0.5,v,'UniformOutput',false)}];
    ub = [{cellfun(@(it)rand(size(it))+0.5,x,'UniformOutput',false)};{cellfun(@(it)rand(size(it))+0.5,v,'UniformOutput',false)}];
    tau = [{cellfun(@(it)rand(size(it)),x,'UniformOutput',false)};{cellfun(@(it)rand(size(it)),v,'UniformOutput',false)}];

else
	lb = {cellfun(@(it)rand(size(it))-0.5,x,'UniformOutput',false)};
    ub = {cellfun(@(it)rand(size(it))+0.5,x,'UniformOutput',false)};
    tau = {cellfun(@(it)rand(size(it)),x,'UniformOutput',false)};
end
    
merit = @(f,dE,bE,varargin) l1merit(f,dE,bE,ub,lb,rho,tau,varargin{:});





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


maxError = max(e);


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

