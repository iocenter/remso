
clc
clear
clear global

% Include REMSO functionalities
addpath(genpath('../../optimization/multipleShooting'));
addpath(genpath('../../optimization/plotUtils'));
addpath(genpath('../../optimization/remso'));
addpath(genpath('../../optimization/remsoSequential'));
addpath(genpath('../../optimization/remsoCrossSequential'));
addpath(genpath('../../optimization/parallel'));
addpath(genpath('../../optimization/testFunctions'));
addpath(genpath('../../optimization/utils'));

addpath(genpath('../../mrstDerivated'));
addpath(genpath('../../mrstLink'));
addpath(genpath('../../mrstLink/wrappers/procedural'));
addpath(genpath('../../netLink'));
addpath(genpath('../../netLink/plottings'));


addpath(genpath('model'));


totalPredictionSteps = 100;
totalControlSteps = 50;

simPerCtrl = totalPredictionSteps/totalControlSteps;
stepControl = cell2mat(arrayfun(@(index)ones(simPerCtrl,1)*index,(1:totalControlSteps)','UniformOutput',false));
ci = @(kk)controlIncidence(stepControl,kk);


W = repmat(struct('cells'   , 0,           ...
    'type'    , 'bhp',             ...
    'val'     , 1*barsa,              ...
    'r'       , 1*meter,           ...
    'dir'     , [],              ...
    'WI'      , 1,                   ...
    'dZ'      , 1, ...
    'name'    , 'Thiago',             ...
    'compi'   , [1 0 0],           ...
    'refDepth', 0,         ...
    'lims'    , [],         ...
    'sign'    , -1,             ...+--i
    'status'  , true,                    ...
    'cstatus' , true),3,1);

W(2).name = 'codas';
W(2).sign = 1;
W(3).name = 'sthener';

%%TODO: pass reservoir state
wellSol = initWellSolLocal(W, reservoirState);

A = [0.1,0.2;
     0.3,0.4];
 
nx = size(A,1);

nw = 2; 
B = sparse(1:nw,1:nw,ones(nw,1),nx,nw);

C = 0.1*ones(3*nw,nx);

D = 0.5*ones(3*nw,nw);

cx = [5,6];
cv = (1:3*nw)/nw;
cu = 1:nw;

Qx = [1,0;0,2];
Qu = eye(nw);
Qv = eye(3*nw);



state = [0.1;0.2];

u1 = ones(nw,1);

netSol = prodNetwork(wellSol);
nScale = [5*barsa;5*barsa];

dpChokes = arroba(@chokesDp,[1,2],{netSol, nScale}, true);

ss.step = repmat({@(xS,u,varargin) linearModel(xS,u,A,B,C,D,'algFun',[],varargin{:})},totalPredictionSteps,1);
ss.ci = ci;
ss.state = state;

bias = 0;
scale = 1;%1e-3;


%obj = @(x,u,v,varargin) linearObjective(x,u,v,cx,cu,cv,varargin{:});
obj = @(x,u,v,varargin) quadraticObjective(x,u,v,cx,cu,cv,Qx,Qu,Qv,'bias',bias,'scale',scale,varargin{:});


lbx = repmat({[-1;-2]},totalPredictionSteps,1);
ubx = repmat({[1;1]},totalPredictionSteps,1);

lbu = repmat({-4},totalControlSteps,1);
ubu = repmat({4},totalControlSteps,1);

lbv = repmat({-10},totalPredictionSteps,1);
ubv = repmat({10},totalPredictionSteps,1);

%{
lbx = repmat({[-inf;-inf]},totalPredictionSteps,1);
ubx = repmat({[inf;inf]},totalPredictionSteps,1);

lbu = repmat({-inf},totalControlSteps,1);
ubu = repmat({inf},totalControlSteps,1);

lbv = repmat({-inf},totalPredictionSteps,1);
ubv = repmat({inf},totalPredictionSteps,1);
%}

u = repmat({u1},totalControlSteps,1);




withAlgs = true;
%[ errorMax,eCross ] = unitTest(u,ss,obj,'totalSteps',[],'debug',true,'feasible',false,'noise',true)





u0 = u;
[~,~,~,simVars,x0,v0] = simulateSystemSS(u0,ss,[]);


[~,~,xd,vd,ax,Ax,av,Av] = condensing(x0,u0,v0,ss,'computeCorrection',true,'withAlgs',withAlgs,'simVars',simVars);

assert(norm(cell2mat(xd),inf)<sqrt(eps))
assert(norm(cell2mat(vd),inf)<sqrt(eps))

%%% if problem is linear, the particular choise of these shouldn't affect
%%% the hessian
muU = cellfun(@(zi)rand(size(zi))',u0','UniformOutput',false);
muX = cellfun(@(zi)rand(size(zi))',x0','UniformOutput',false);
muV = cellfun(@(zi)rand(size(zi))',v0','UniformOutput',false);


[f0,objPartials] = obj(x0,u0,v0,'gradients',true);
% gbar = g+nu
gbar.Jx =  cellfun(@(Jz,mi)(Jz+mi),objPartials.Jx,muX,'UniformOutput',false);
gbar.Ju =  cellfun(@(Jz,mi)(Jz+mi),objPartials.Ju,muU,'UniformOutput',false);
if withAlgs
    gbar.Jv = cellfun(@(Jz,mi)(Jz+mi),objPartials.Jv,muV,'UniformOutput',false);
end
[~,~,~,lambdaX,lambdaV]= simulateSystemZ(u0,x0,v0,ss,[],'simVars',simVars,'JacTar',gbar,'withAlgs',withAlgs);


%[ hhxvu ] = buildFullHessian( x0,v0,u0,obj,ss,lambdaX,lambdaV,'withAlgs',withAlgs);

uDims = cellfun(@(ui)numel(ui),u0);
xDims = cellfun(@(ui)numel(ui),x0);
vDims = cellfun(@(ui)numel(ui),v0);



%nullspace and range space

Z = [cell2matFill(Ax);cell2matFill(Av);eye(sum(uDims))];
Y = [speye(sum(xDims+vDims));zeros(sum(uDims),sum([xDims;vDims]))];


QX = arrayfun(@(x)Qx,xDims,'UniformOutput',false);
QU = arrayfun(@(x)Qu,uDims,'UniformOutput',false);
QV = arrayfun(@(x)Qv,vDims,'UniformOutput',false);

Q = blkdiag(QX{:},QV{:},QU{:});


%norm(cell2mat(hhxvu)-Q,inf)  %% just to check

B = cell2mat([arrayfun(@(x)cx,xDims','UniformOutput',false),...
              arrayfun(@(x)cv,vDims','UniformOutput',false),...
              arrayfun(@(x)cu,uDims','UniformOutput',false)]);
          
[xs,vs,Jac,convergence,simVars,usliced] = simulateSystem(x0,u0,ss,'gradients',true,'withAlgs',withAlgs);

          
xJx = cell2matFill(Jac.xJx,xDims,xDims')-speye(sum(xDims));
xJu = cell2matFill(Jac.xJu,xDims,uDims');
vJx = cell2matFill(Jac.vJx,vDims,xDims');
vJu = cell2matFill(Jac.vJu,vDims,uDims');

xJv = sparse(sum(xDims),sum(vDims));
vJv = -speye(sum(vDims));


A = [xJx,xJv,xJu;
     vJx,vJv,vJu];

z0 = cell2mat([x0;v0;u0]);

b = A*z0;


dB = B+ z0'*Q;


LB = cell2mat([lbx;lbv;lbu])-z0;
UB = cell2mat([ubx;ubv;ubu])-z0;


[dz,~,~,~,lambda] = quadprog(Q,dB,[],[],A,zeros(size(A,1),1),LB,UB);

z = z0 + dz;


z = mat2cell(z,[xDims;vDims;uDims],1);
xOpt = z(1:numel(x0));
vOpt = z(numel(x0)+1:numel(x0)+numel(v0));
uOpt = z(numel(x0)+numel(v0)+1:end);



[f,~,~,simVars,x0,v0] = simulateSystemSS(uOpt,ss,obj);
norm(cell2mat(xOpt)-cell2mat(x0),inf)
norm(cell2mat(vOpt)-cell2mat(v0),inf)


bias = -f/scale;
obj = @(x,u,v,varargin) quadraticObjective(x,u,v,cx,cu,cv,Qx,Qu,Qv,'bias',bias,'scale',scale,varargin{:});



Mopt = scale*Z'*Q*Z;

u = cellfun(@(zi)rand(size(zi)),u0,'UniformOutput',false);
x = cellfun(@(zi)rand(size(zi)),x0,'UniformOutput',false);
v = cellfun(@(zi)rand(size(zi)),v0,'UniformOutput',false);


[u,x,v,f,xd,M,simVars] = remso(u,ss,obj,'lbx',lbx,'ubx',ubx,'lbu',lbu,'ubu',ubu,'lbv',lbv,'ubv',ubv,...
    'lkMax',4,'max_iter',200,'tol',1e-5,'plot',false,'x',x,'v',v,'condense',false,'M',Mopt,...
    'skipRelaxRatio',inf);




u1 = cellfun(@(zi)rand(size(zi)),u0,'UniformOutput',false);
x1 = cellfun(@(zi)rand(size(zi)),x0,'UniformOutput',false);
v1 = cellfun(@(zi)rand(size(zi)),v0,'UniformOutput',false);

%[f,~,~,simVars,xs,vs] = simulateSystemSS(u,ss,obj);
[f0,objPartials1] = obj(x1,u1,v1,'gradients',true);


u2 = cellfun(@(zi)rand(size(zi)),u0,'UniformOutput',false);
x2 = cellfun(@(zi)rand(size(zi)),x0,'UniformOutput',false);
v2 = cellfun(@(zi)rand(size(zi)),v0,'UniformOutput',false);

%[f,~,~,simVars,xs,vs] = simulateSystemSS(u,ss,obj);
[f0,objPartials2] = obj(x2,u2,v2,'gradients',true);


muU = cellfun(@(zi)rand(size(zi))',u0','UniformOutput',false);
muX = cellfun(@(zi)rand(size(zi))',x0','UniformOutput',false);
muV = cellfun(@(zi)rand(size(zi))',v0','UniformOutput',false);


gbar1.Jx =  cellfun(@(Jz,mi)(Jz+mi),objPartials1.Jx,muX,'UniformOutput',false);
gbar1.Ju =  cellfun(@(Jz,mi)(Jz+mi),objPartials1.Ju,muU,'UniformOutput',false);
gbar1.Jv = cellfun(@(Jz,mi)(Jz+mi),objPartials1.Jv,muV,'UniformOutput',false);

gbar2.Jx =  cellfun(@(Jz,mi)(Jz+mi),objPartials2.Jx,muX,'UniformOutput',false);
gbar2.Ju =  cellfun(@(Jz,mi)(Jz+mi),objPartials2.Ju,muU,'UniformOutput',false);
gbar2.Jv = cellfun(@(Jz,mi)(Jz+mi),objPartials2.Jv,muV,'UniformOutput',false);


[~,gbarZ1] = simulateSystemZ(u1,x1,v1,ss,[],'JacTar',gbar1);
[~,gbarZ2] = simulateSystemZ(u2,x2,v2,ss,[],'JacTar',gbar2);


[~,gZ1] = simulateSystemZ(u1,x1,v1,ss,[],'JacTar',objPartials1);
[~,gZ2] = simulateSystemZ(u2,x2,v2,ss,[],'JacTar',objPartials2);

assert(norm(cell2mat(gbarZ1)-cell2mat(gbarZ2) - (cell2mat(gZ1) - cell2mat(gZ2) ),inf) < sqrt(eps))  %% the lagrangian is independent of mu




uBV = cell2mat(u);
uV = cell2mat(u);

nru = numel(uV);
M = eye(nru);

[f,gbarZm,~,simVars,x0,v0] = simulateSystemSS(mat2cell(uV,uDims,1),ss,obj,'gradients',true);

S= [];
Y =[];
it = 1;
[U,S,V] = svd(M-Mopt);
sMax = S(1,1);

while sMax > sqrt(eps) && it < 1000
    
    uV = uBV+U(:,1); 
    %ek = zeros(size(uBV));
    %ek(it) = ek(it) +1;
    %ek = -M\cell2mat(gbarZm)';
    %uV = uBV+ek;
    
    [f,gbarZ,~,simVars,x0,v0] = simulateSystemSS(mat2cell(uV,uDims,1),ss,obj,'gradients',true);

	y = cellfun(@(gbarZi,gbarZmi)gbarZi-gbarZmi,gbarZ,gbarZm,'UniformOutput',false);
	y = cell2mat(y);
	s = uV-uBV;
        
    [ M,skipping,damping,minEig,sTy,theta] = dampedBfgsUpdate(M,y,s,'allowDamp',false);    

    gbarZm = gbarZ;
    uBV = uV;
    
    [U,S,V] = svd(M-Mopt);
    sMax = S(1,1);
    [sTy,sMax]
    
    it = it +1;
    %diag(S)
end

