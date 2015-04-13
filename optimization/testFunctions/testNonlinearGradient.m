function [ maxError,eCross ] = testNonlinearGradient(x,u,v,ss,obj,varargin)


opt = struct('pert',1e-5,'debug',false,'withAlgs',false);
opt = merge_options(opt, varargin{:});

maxError = -inf;
pertV = opt.pert;

xDims = cellfun(@(z)numel(z),x);
uDims = cellfun(@(z)numel(z),u);
vDims = cellfun(@(z)numel(z),v);

withAlgs = opt.withAlgs;


%x = cellfun(@(zi)zi+norm(zi)*rand(size(zi))/20,x,'UniformOutput',false);
%v = cellfun(@(zi)zi+norm(zi)*rand(size(zi))/20,v,'UniformOutput',false);
%u = cellfun(@(zi)zi+norm(zi)*rand(size(zi))/20,u,'UniformOutput',false);



[xs,vs,xd,vd,a,A,av,Av] = condensing(x,u,v,ss,'computeCorrection',true,'withAlgs',withAlgs);


[f,gradU] = targetGrad(x,u,v,obj,A,Av,ss.ci);


%{
%% this is being tested elsewhere

[xsF,vF,JacFull] = simulateSystem(x,u,ss,'gradients',true,'withAlgs',withAlgs);

xJu = cell2matFill(JacFull.xJu,[nx,nu]);
xJx = cell2matFill(JacFull.xJx,[nx,nx]);
vJu = cell2matFill(JacFull.vJu,[nv,nu]);
vJx = cell2matFill(JacFull.vJx,[nv,nx]);




AB = [xJx-eye(size(xJx)),zeros(size(xJx,1),size(vJx,1)),xJu;
      vJx,-eye(size(vJx,1)),vJu                                    ];


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




[fz,Bz,simVars,usz,lambdaX,lambdaV] = simulateSystemZ(u,x,v,ss,obj,'withAlgs',withAlgs);






e = [];
e = [e f-fz norm(cell2mat(gradU)-cell2mat(Bz))];



%test JacTar in simulateSystemZ 
testObj.Jx = cellfun(@(zi)rand([size(zi,2),size(zi,1)]),x','UniformOutput',false);
testObj.Jv = cellfun(@(zi)rand([size(zi,2),size(zi,1)]),v','UniformOutput',false);
testObj.Ju = cellfun(@(zi)rand([size(zi,2),size(zi,1)]),u','UniformOutput',false);


if withAlgs
    gZ = vectorTimesZ(testObj.Jx,testObj.Ju,testObj.Jv,A,Av,ss.ci );
else
    gZ = vectorTimesZ(testObj.Jx,testObj.Ju,[],A,[],ss.ci );
end


[xsR,vsR,~,convergedR,simVarsR,uslicedR] = simulateSystem(x,u,ss,'gradients',false,'guessX',xs,'guessV',vs,'withAlgs',withAlgs);



[~,gradUY,~,~,~,~] = simulateSystemZ(u,x,v,ss,obj,'simVars',simVarsR,'JacTar',testObj,'withAlgs',withAlgs);



wErr = cellfun(@minus,gradUY,gZ,'UniformOutput',false);


e = [e norm(cell2mat(wErr))];





activeSet.lb.u = arrayfun(@(udi)rand(udi,1)<0.5,uDims,'UniformOutput',false);
activeSet.ub.u = arrayfun(@(udi)rand(udi,1)<0.5,uDims,'UniformOutput',false);

activeSet.lb.x = arrayfun(@(udi)rand(udi,1)<0.5,xDims,'UniformOutput',false);
activeSet.ub.x = arrayfun(@(udi)rand(udi,1)<0.5,xDims,'UniformOutput',false);

activeSet.lb.v = arrayfun(@(udi)rand(udi,1)<0.5,vDims,'UniformOutput',false);
activeSet.ub.v = arrayfun(@(udi)rand(udi,1)<0.5,vDims,'UniformOutput',false);


[ t,Jac ] = activeSet2TargetXV(uDims,activeSet );


[~,Aact,~,~,~,~] = simulateSystemZ(u,x,v,ss,obj,'simVars',simVarsR,'JacTar',Jac,'withAlgs',withAlgs);


Axm = cell2matFill(A,xDims,uDims);
Avm = cell2matFill(Av,vDims,uDims);

AmAct = [-Axm(cell2mat(activeSet.lb.x),:);
              Axm(cell2mat(activeSet.ub.x),:);
             -Avm(cell2mat(activeSet.lb.v),:);
              Avm(cell2mat(activeSet.ub.v),:)];
Aact = cell2mat(Aact);
          
eM = AmAct-Aact;

e = [e norm(eM,inf)];


%% test LagrangianF


% this affect the cross-term approx more than expected
lambdaX = cellfun(@(zi)rand(size(zi))',x','UniformOutput',false);
lambdaV = cellfun(@(zi)rand(size(zi))',v','UniformOutput',false);

muU = cellfun(@(zi)rand(size(zi))',u','UniformOutput',false);
muX = cellfun(@(zi)rand(size(zi))',x','UniformOutput',false);
muV = cellfun(@(zi)rand(size(zi))',v','UniformOutput',false);

lbx = cellfun(@(zi)min(zi+rand(size(zi))-0.5,zi),x,'UniformOutput',false);
ubx = cellfun(@(zi,lbxi)max(max(zi+rand(size(zi))-0.5,lbxi+0.1),zi),x,lbx,'UniformOutput',false);

lbv = cellfun(@(zi)min(zi+rand(size(zi))-0.5,zi),v,'UniformOutput',false);
ubv = cellfun(@(zi,lbxi)max(max(zi+rand(size(zi))-0.5,lbxi+0.1),zi),v,lbv,'UniformOutput',false);

lbu = cellfun(@(zi)min(zi+rand(size(zi))-0.5,zi),u,'UniformOutput',false);
ubu = cellfun(@(zi,lbxi)max(max(zi+rand(size(zi))-0.5,lbxi+0.1),zi),u,lbu,'UniformOutput',false);



[ lagF,lagG] = lagrangianF( u,x,v,lambdaX,lambdaV,muU,muX,muV,obj,ss,lbx,lbv,lbu,ubx,ubv,ubu,'gradients',true,'withAlgs',withAlgs);

% test lagG
lagFun = @(dd) lagrangianF(mat2cell(dd,uDims,1),x,v,lambdaX,lambdaV,muU,muX,muV,obj,ss,lbx,lbv,lbu,ubx,ubv,ubu,'gradients',false);
lagJu = calcPertGrad(lagFun,cell2mat(u),opt.pert);

lagFun = @(dd) lagrangianF(u,mat2cell(dd,xDims,1),v,lambdaX,lambdaV,muU,muX,muV,obj,ss,lbx,lbv,lbu,ubx,ubv,ubu,'gradients',false);
lagJx = calcPertGrad(lagFun,cell2mat(x),opt.pert);

lagFun = @(dd) lagrangianF(u,x,mat2cell(dd,vDims,1),lambdaX,lambdaV,muU,muX,muV,obj,ss,lbx,lbv,lbu,ubx,ubv,ubu,'gradients',false);
lagJv = calcPertGrad(lagFun,cell2mat(v),opt.pert);

% test lagrangian function
 e = [e, max(abs(cell2mat(lagG.Ju)'-lagJu')),max(abs(cell2mat(lagG.Jx)'-lagJx')),max(abs(cell2mat(lagG.Jv)'-lagJv'))];

 
 [lagG2] = lagrangianG(  u,x,v,lambdaX,lambdaV,muU,muX,muV,obj,ss,'withAlgs',withAlgs);
 
 e = [e, max(abs(cell2mat(lagG.Ju)-cell2mat(lagG2.Ju))),max(abs(cell2mat(lagG.Jx)-cell2mat(lagG2.Jx))),max(abs(cell2mat(lagG.Jv)-cell2mat(lagG2.Jv)))];

crossAlgo = nargin(@l1merit)  ~= -8 ;

eCross  = testCrossTerm( x,v,u,obj,ss,muX,muV,muU,withAlgs );




uSeed = cellfun(@(it)rand(size(it)),u,'UniformOutput',false);

[~,~,~,~,~,AxS,~,AvS] = condensing(x,u,v,ss,'uRightSeeds',uSeed,'computeNullSpace',true,'withAlgs',withAlgs);

e = [e, norm(cell2matFill(A,xDims,uDims')*cell2mat(uSeed)-cell2mat(AxS)),norm(cell2matFill(Av,vDims,uDims')*cell2mat(uSeed)-cell2mat(AvS)) ];



xdRand = cellfun(@(zi)rand(size(zi)),x,'UniformOutput',false);
vdRand = cellfun(@(zi)rand(size(zi)),v,'UniformOutput',false);


[~,~,~,~,axRand,~,avRand,~] = condensing(x,u,v,ss,'computeCorrection',true,'computeNullSpace',false,'xd',xdRand,'vd',vdRand,'withAlgs',withAlgs);

auRand = cellfun(@(ui)zeros(size(ui)),u,'UniformOutput',false);

[~,~,JacRand,~,~,~] = simulateSystem(x,u,ss,'gradients',true,'guessX',xs,'guessV',vs,'xRightSeed',axRand,'uRightSeed',auRand,'withAlgs',withAlgs);

e = [e,norm(cell2mat(xdRand) - (cell2mat(axRand) - cell2mat(JacRand.xJ)))];
e = [e,norm(cell2mat(vdRand) - (cell2mat(avRand) - cell2mat(JacRand.vJ)))];


[f,JacTarFull] = obj(xs,u,vs,'gradients',true);


fx = @(xit) obj(mat2cell(xit,xDims,1),u,vs);
dfdx = calcPertGrad(fx,cell2mat(xs),pertV);


fu = @(uit) obj(xs,mat2cell(uit,uDims,1),vs);
dfdu = calcPertGrad(fu,cell2mat(u),pertV);

fv = @(vit) obj(xs,u,mat2cell(vit,vDims,1));
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
    gradUF = cell2mat(JacTarFull.Ju) + cell2mat(JacTarFull.Jx) * cell2matFill(A,xDims,uDims') + cell2mat(JacTarFull.Jv) * cell2matFill(Av,vDims,uDims');
else
    gradUF = cell2mat(JacTarFull.Ju) + cell2mat(JacTarFull.Jx) * cell2matFill(A,xDims,uDims') ;
end
e = [e norm(cell2mat(gradU)-gradUF)];

if ~crossAlgo
    [ em ] = testMerit( x,v,u,obj,ss,withAlgs);
else
    [ em ] = testMeritCross( x,v,u,obj,ss,withAlgs);
end

e = [e em];

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

