function [ maxerror ] = testLiftOptReduction_R(x,u,v,ss)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

nR = numel(x);


[xs,vs,s2,JacMS,converged,simVars,usliced] = simulateSystem_R(x,u,v,ss,'gradients',true);

s = rand(numel(s2),1);

xDims = cellfun(@(z)cellfun(@numel,z),x,'UniformOutput',false);
uDims = cellfun(@(z)numel(z),u);
vDims = cellfun(@(z)cellfun(@numel,z),v,'UniformOutput',false);


xd = cell2mat(cellfun(@(xsFr,xr)cell2mat(cellfun(@minus,xsFr,xr,'UniformOutput',false)),xs,x,'UniformOutput',false));
vd = cell2mat(cellfun(@(xsFr,xr)cell2mat(cellfun(@minus,xsFr,xr,'UniformOutput',false)),vs,v,'UniformOutput',false));
sd = s2-s;
xs = cell2mat(cellfun(@cell2mat,xs,'UniformOutput',false));
vs = cell2mat(cellfun(@cell2mat,vs,'UniformOutput',false));


xJx = JacMS.xJx;
xJx = cellfun(@(xdJacr,xDimsr)cell2matFill(xdJacr,xDimsr,xDimsr),xJx,xDims,'UniformOutput',false);
xJx = blkdiag(xJx{:});
xdJac = xJx - speye(size(xJx,1));

xJu = cellfun(@(xdJacr,xDimsr)cell2matFill(xdJacr,xDimsr,uDims),JacMS.xJu,xDims,'UniformOutput',false);
xJu = cell2mat(xJu);

vJx = JacMS.vJx;
vJx = cellfun(@(vJxr,vDimsr,xDimsr)cell2matFill(vJxr,vDimsr,xDimsr),vJx,vDims,xDims,'UniformOutput',false);
vJx = blkdiag(vJx{:});

vJu = JacMS.vJu;
vJu = cellfun(@(vJxr,vDimsr)cell2matFill(vJxr,vDimsr,uDims),vJu,vDims,'UniformOutput',false);
vJu = cell2mat(vJu);


sJv = JacMS.sJv;
sJx = JacMS.sJx;
sJu = JacMS.sJu;

sJv = cell2mat(cellfun(@cell2mat,sJv,'UniformOutput',false));
sJx = cell2mat(cellfun(@cell2mat,sJx,'UniformOutput',false));
sJu = cell2mat(sJu);


% lift-opt related matrix brute force calculated
invdhmidx = inv(xdJac);
ax = -invdhmidx*xd;
Ax  = -invdhmidx*xJu;

av = (vd + vJx*ax);
Av = (vJu+vJx*Ax);

as = (sd + sJv*av+ sJx*ax);
As = (sJu+sJv*Av+sJx*Ax);

[xsZ,vsZ,s2Z,xdZ,vdZ,sdZ,axZ,AxZ,avZ,AvZ,asZ,AsZ] = condensing_R(x,u,v,s,ss,'computeCorrection',true);

xsZ = cell2mat(cellfun(@cell2mat,xsZ,'UniformOutput',false));
vsZ = cell2mat(cellfun(@cell2mat,vsZ,'UniformOutput',false));

xdZ = cell2mat(cellfun(@cell2mat,xdZ,'UniformOutput',false));
vdZ = cell2mat(cellfun(@cell2mat,vdZ,'UniformOutput',false));

axZ = cell2mat(cellfun(@cell2mat,axZ,'UniformOutput',false));
avZ = cell2mat(cellfun(@cell2mat,avZ,'UniformOutput',false));

AxZ = cellfun(@(xdJacr,xDimsr)cell2matFill(xdJacr,xDimsr,uDims),AxZ,xDims,'UniformOutput',false);
AxZ = cell2mat(AxZ);
AvZ = cellfun(@(xdJacr,xDimsr)cell2matFill(xdJacr,xDimsr,uDims),AvZ,vDims,'UniformOutput',false);
AvZ = cell2mat(AvZ);
AsZ = cell2mat(AsZ);

e1 = norm(xsZ-xs);
e2 = norm(vsZ-vs);
e3 = norm(s2Z-s2);
e4 = norm(xdZ-xd);
e5 = norm(vdZ-vd);
e6 = norm(sdZ-sd);
e7 = norm(axZ-ax);
e8 = norm(AxZ-Ax,inf);
e9 = norm(avZ - av);
e10 = norm(AvZ - Av,inf);
e11 = norm(asZ - as);
e12 = norm(AsZ - As,inf);

xF = cell2mat(cellfun(@cell2mat,x,'UniformOutput',false));
vF = cell2mat(cellfun(@cell2mat,v,'UniformOutput',false));

e13 = norm(xs-(xF+xd));
e14 = norm(vs-(vF+vd));
e15 = norm(s2-(s +sd));

xT = sum(cell2mat(xDims));
vT = sum(cell2mat(vDims));
sT = numel(s);
uT = sum(uDims);

AT = [xdJac,            sparse(xT,vT),  sparse(xT,sT),  xJu;
      vJx,              -speye(vT),     sparse(vT,sT),  vJu;
      sJx,              sJv,            -speye(sT),     sJu];

nullspaceX = [Ax;Av;As;speye(uT)];
    
e16 = norm(AT*nullspaceX,inf);

e17 = norm([xd;vd;sd] + AT*[ax;av;as;zeros(uT,1)]);

nSeeds = 1 + floor(5*rand);

xLeftSeed = cellfun(@(xDimsr)arrayfun(@(xdirk)rand(nSeeds,xdirk),xDimsr','UniformOutput',false),xDims','UniformOutput',false);
vLeftSeed = cellfun(@(xDimsr)arrayfun(@(xdirk)rand(nSeeds,xdirk),xDimsr','UniformOutput',false),vDims','UniformOutput',false);
sLeftSeed = rand(nSeeds,numel(s));


LS = [cell2mat(cellfun(@cell2mat,xLeftSeed,'UniformOutput',false)),...
      cell2mat(cellfun(@cell2mat,vLeftSeed,'UniformOutput',false)),...
      sLeftSeed];

xRightSeed = cellfun(@(xDimsr)arrayfun(@(xdirk)rand(xdirk,nSeeds),xDimsr,'UniformOutput',false),xDims,'UniformOutput',false);
uRightSeed = arrayfun(@(xdirk)rand(xdirk,nSeeds),uDims,'UniformOutput',false);
vRightSeed = cellfun(@(xDimsr)arrayfun(@(xdirk)rand(xdirk,nSeeds),xDimsr,'UniformOutput',false),vDims,'UniformOutput',false);

RS = [cell2mat(cellfun(@cell2mat,xRightSeed,'UniformOutput',false));...
      cell2mat(cellfun(@cell2mat,vRightSeed,'UniformOutput',false));...
      cell2mat(uRightSeed)];


AMS = [xJx,            sparse(xT,vT),	xJu;
       vJx,            sparse(vT,vT),	vJu;
       sJx,            sJv,             sJu];  
  

[xs,vs,s2,JacMS,converged,simVars,usliced] = simulateSystem_R(x,u,v,ss,'gradients',true,'xLeftSeed',xLeftSeed,'vLeftSeed',vLeftSeed,'sLeftSeed',sLeftSeed);


LSAMS = [cell2mat(cellfun(@cell2mat,JacMS.Jx,'UniformOutput',false)),...
         cell2mat(cellfun(@cell2mat,JacMS.Jv,'UniformOutput',false)),...
         cell2mat(JacMS.Ju)];



e18 = norm(LSAMS - LS*AMS,inf);


[xs,vs,s2,JacMS,converged,simVars,usliced] = simulateSystem_R(x,u,v,ss,'gradients',true,'xRightSeed',xRightSeed,'uRightSeed',uRightSeed,'vRightSeed',vRightSeed);


AMSRS = [cell2mat(cellfun(@cell2mat,JacMS.xJ,'UniformOutput',false));...
         cell2mat(cellfun(@cell2mat,JacMS.vJ,'UniformOutput',false));...
         JacMS.sJ];

e19 = norm(AMSRS -AMS*RS,inf);




[obj,JacTar] = targetAll(x,u,v,s);


[gU,lambdaX,lambdaV,lambdaS] = simulateSystemZ_R(u,x,v,ss,JacTar);

lambdaX = cell2mat(cellfun(@cell2mat,lambdaX,'UniformOutput',false));
lambdaV = cell2mat(cellfun(@cell2mat,lambdaV,'UniformOutput',false));


gUc = cell2mat(cellfun(@cell2mat,JacTar.Jx,'UniformOutput',false))*Ax + ...
      cell2mat(cellfun(@cell2mat,JacTar.Jv,'UniformOutput',false))*Av + ...
      JacTar.Js*As +...
      cell2mat(JacTar.Ju);
  
e19 = norm(cell2mat(gU)-gUc,'inf');



Y = [speye(xT+vT+sT);sparse(uT,xT+vT+sT)];

gx = cell2mat(cellfun(@cell2mat,JacTar.Jx,'UniformOutput',false));
gv = cell2mat(cellfun(@cell2mat,JacTar.Jv,'UniformOutput',false));
gs = JacTar.Js;
gu = cell2mat(JacTar.Ju);


e20 = norm((-(Y'*AT')\(Y'*[gx';gv';gs';gu'])) - [lambdaX,lambdaV,lambdaS]');




maxerror = max([e1,e2,e3,e4,e5,e6,e7,e8,e9,e10,e11,e12,e13,e14,e15,e16,e17,e18,e19]);

end

