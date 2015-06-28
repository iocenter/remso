function [u,x,v,s,f,M,simVars] = remso(u,sss,obj,varargin)
% REMSO
% REservoir Multiple Shooting Optimization.
% REduced Multiple Shooting Optimization.
%
% This is the main interface to the REMSO solver.
%
% SYNOPSIS:
%  [u,x,v,f,xd,M,simVars] = remso(u, ss, obj)
%  [u,x,v,f,xd,M,simVars] = remso(u, ss, obj, 'pn', pv, ...)
% PARAMETERS:
%   u - cellarray containing a initial control guess for each control
%       period.
%
%   ss - A simulator structure, containing all the required
%        information on the model.
%
%   obj - A nonlinear function structure defining the objective function
%
%   'pn'/pv - List of 'key'/value pairs defining optional parameters. The
%             supported options are:
%
%   lbx - State lower bound for each point in the prediction horizon.
%
%   ubx - State upper bound for each point in the prediction horizon.
%
%   lbv - Algebraic state lower bound for each point in the prediction horizon.
%
%   ubv - Algebraic state upper bound for each point in the prediction horizon.
%
%   lbu - Control input lower bound for each control period.
%
%   ubu - Control input upper bound for each control period.
%
%   tol - Master tolerance.
%
%   tolU - Convergence tolerance for the controls.
%
%   tolX - Convergence tolerance for the states.
%
%   tolV - Convergence tolerance for the algebraic variables.
%
%   max_iter - Maximum iterations allowed for the main algorithm.
%
%   M - Initial reduced hessian approximation.
%
%   x - Initial guess for the states in the prediction horizon..
%
%   v - Initial guess for the algebraic states in the control horizon.
%
%   plotFunc - plotFunc(x,u,v,xd).  Plot function for the current solution
%              iterate.
%
%   lkMax - Maximum number of evaluated points during line-search.
%
%   eta - Constant related to the Wolf curvature condition.
%
%   tauL - Constant related to the minimum descent condition.
%
%   debugLS - Plot debug information during line-search.
%
%   qpDebug - Print debug information related to the QP solving process.
%
%   lowActive - Initial active set estimate related to the lower bounds.
%
%   upActive - Initial active set estimate related to the upper bounds.
%
%   simVars - Simulation variables, for hot start initialization.
%
%   debug - Print debug information containing general algorithm
%           performance.
%
%   plot - Flag to allow plotting at each iteration.
%
%   saveIt - Save current iterate variables at each iteratoin.
%
%
% RETURNS:
%
%   u - Optimal control estimate.
%
%   x - State forecast.
%
%   v - Algebraic state forecast.
%
%   f - Estimated objective function value.
%
%   xd - State forecast error estimation.
%
%   M - Hessian approximation.
%
%   simVars - Final simulation variables.
%
% SEE ALSO:
%
%
%{

Copyright 2013-2014, Andres Codas.

REMSO is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

REMSO is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with REMSO.  If not, see <http://www.gnu.org/licenses/>.

%}
opt = struct('lbx',[],'ubx',[],'lbv',[],'ubv',[],'lbu',[],'ubu',[],...
    'ubs',[],'lbs',[],...
    'tol',1e-1,'tolU',1e-2,'tolX',1e-2,'tolV',1e-2,'tolS',1e-2,'max_iter',50,...
    'M',[],'x',[],'v',[],...
    'plotFunc',[],...
    'BFGSRestartscale', true,'BFGSRestartmemory',6,...
    'lkMax',4,'eta',0.1,'tauL',0.1,'debugLS',false,'curvLS',true,...
    'rhoHat',1e-10,...
    'nCons',0.25,...  % fraction of constraints to be added to the QP problem, if inf intended, then set it to 0;
    'qpDebug',true,...
    'lowActive',[],'upActive',[],...
    'simVars',[],'debug',true,'plot',false,'saveIt',false,...
    'controlWriter',[],...
    'multiplierFree',inf,...
    'allowDamp',true,...
    'qpFeasTol',1e-6,...
    'computeCrossTerm',true,...
    'condense',false,'testQP',false);

opt = merge_options(opt, varargin{:});


masterTol = min([opt.tol,opt.tolU,opt.tolX,opt.tolV]);

%The qpFeasTol must be tighter than tol, tolX, tolV, and tolU'
if opt.qpFeasTol > masterTol
    opt.qpFeasTol = masterTol;
end

debug = opt.debug;

% number of variables
ss = sss.ss;
if opt.debug
	%spmd
		outputName = sprintf('w%d.log', sss.jobSchedule.my_rank+1);
		fidW = fopen(outputName,'w');
	%end
else
	fidW = [];
end
sss.jobSchedule.fidW = fidW;
jobSchedule = sss.jobSchedule;
imMaster = jobSchedule.imMaster;


%spmd
xDims0 = getXDims0(ss);  %% it is assumed that the realizations can have different number of states
%end
uDims = cellfun(@numel,u);

% dimension of the control space, dimension of the reduced problem
nru = sum(uDims);
if ~isfield(sss,'eta')
    sss.eta = 0;
    warning('sss.eta not specified. Adopted sss.eta = 0, i.e., mean value as risk measure')
end

nCons = ceil(nru*opt.nCons);
%% Control and state bounds processing
lbu = opt.lbu;
if isempty(lbu)
    lbu = buildRelaxedBoundS(uDims,-1);
end
ubu = opt.ubu;
if isempty(ubu)
    ubu = buildRelaxedBoundS(uDims,1);
end

[~,u]  = checkBounds(     lbu,u,    ubu,'chopp',true,'verbose',debug);
uV = cell2mat(u);

% contrasting the other implementations, use a single lbx for all the
% prediction horizon.
lbx = opt.lbx;
if isempty(lbx)
    %spmd
    lbx =  buildRelaxedBoundS(xDims0,-1);
    %end
end
ubx = opt.ubx;
if isempty(ubx)
    %spmd
    ubx =  buildRelaxedBoundS(xDims0,1);
    %end
end
%% initial simulation profile
simVars = opt.simVars;
if isempty(simVars)
    simVars = createEmptyDistributedVar(jobSchedule);  % each cell inside most
end

%% Process initial MS simulation guess, if not given, get it by forward simulation
simulateSS = false;
if ~isempty(opt.x)
    %  Initial guess for prediction given by the user
	x = opt.x;
	[x] = choppBounds( lbx,x,ubx,debug);
	xs = x;
else
    % Initial guess not provided, take from a simulation in the gradient
    % routine
    simulateSS = true;
    x = createEmptyDistributedVar(jobSchedule);
    xs = createEmptyDistributedVar(jobSchedule);
end

if isempty(opt.v)
    v =  createEmptyDistributedVar(jobSchedule);
    vs = createEmptyDistributedVar(jobSchedule);
else
    vs = opt.v;
    v  = opt.v;
	if ~isempty(opt.lbv) || ~isempty(opt.lbv)
        lbv=opt.lbv;
        ubv=opt.ubv;
        %spmd
    	[v] = choppBounds(lbv,v,ubv,debug);
        %end
	end
end

if simulateSS
    [~,~,~,simVars,xs,vs,s2,usliced] = simulateSystemSS_R(u,sss,[],'guessX',xs,'guessV',vs,'simVars',simVars);
    %spmd
    [x,ok]=choppBounds(lbx,xs,ubx,debug);
    ok = gopMPI('*',ok+0,jobSchedule)==1;
    %end
    v = vs;
	okv = true;
	if ~isempty(opt.lbv) || ~isempty(opt.lbv)
        lbv=opt.lbv;
        ubv=opt.ubv;
        %spmd
    	[v,okv] = choppBounds( lbv,v,ubv,debug);
        okv = gopMPI('*',okv+0,jobSchedule)==1;
        %end
	end
	ok = ok && okv;
	ok = NMPI_Bcast(ok+0,1,jobSchedule.Master_rank,jobSchedule.my_rank)==1;
    if ~ok  %% simVars is not correct!
    	[xs,vs,s2,~,~,simVars,usliced] = simulateSystem_R(x,u,v,sss,'gradients',false,'guessX',xs,'guessV',vs);
    end
else
    [xs,vs,s2,~,~,simVars,usliced] = simulateSystem_R(x,u,v,sss,'gradients',false,'guessX',xs,'guessV',vs,'simVars',simVars);
end
s = s2;

%spmd
xDims = getZDims(x);
vDims = getZDims(v);
%end
if imMaster
sDims = numel(s);
end
%lets assume the user is smart enough
%assert(sum([0;cell2mat(vDims)])>0,'The robust optimization formulation must contain algebraic variables to be well defined');

lbv = opt.lbv;
if isempty(lbv)
    %spmd
    lbv = buildRelaxedBound(vDims,-1);
    %end
end
ubv = opt.ubv;
if isempty(ubv)
    %spmd
    ubv = buildRelaxedBound(vDims,1);
    %end
end
ubs = opt.ubs;
if imMaster && isempty(ubs)
    ubs = inf(size(s));
end
lbs = opt.lbs;
if imMaster && isempty(lbs)
    lbs = -inf(size(s));
end

if imMaster
[~,s]=checkBounds(lbs,s,ubs,'chopp',true,'verbose',debug);
end

Aact1= [];
predictor = [];
constraintBuilder = [];

% Multiple shooting simulation function
simFunc = @(xk,uk,vk,varargin) simulateSystem_R(xk,uk,vk,sss,varargin{:});


%% Define empty active sets if they are not given
lowActive = opt.lowActive;
if isempty(lowActive)
    %spmd
    alx = startEmptyActiveSet(xDims);
    alv = startEmptyActiveSet(vDims);
    %end
    lowActive.x = alx;
    lowActive.v = alv;
    if imMaster
    als = {false(sDims,1)};
    lowActive.s = als;
    end
end
upActive = opt.upActive;
if isempty(upActive)
    %spmd
    aux = startEmptyActiveSet(xDims);
    auv = startEmptyActiveSet(vDims);
    %end
    upActive.x = aux;
    upActive.v = auv;
    if imMaster
    aus = {false(sDims,1)};
    upActive.s = aus;
    end
end



%% lagrange multipliers estimate initilization
%spmd
mudx=  initDualVariable(xDims);
mudv = initDualVariable(vDims);
%end
mudu = initDualVariableS(uDims);
if imMaster
muds = initDualVariableS(sDims);
else
muds = nan;
end
%% Hessian Initializaiton
if(isempty(opt.M))
    hInit = true;
    M = eye(nru);
else
    hInit = false;
    M = opt.M;
end

% clean debug file
if debug
    fid = fopen('logBFGS.txt','w');
    fclose(fid);
end


%% Curvature history record
y = zeros(1,sum(uDims));
su = zeros(sum(uDims),1);
sTy = 0;

S = [];
Y = [];



%% Line-search parameters
rho = opt.rhoHat;
rhoHat = opt.rhoHat;
returnVars = [];
relax = false;   % to avoid the hessian update and perform a fine line-search
errorSumB = [];
dualApproxB = [];



% convergence flag
converged = false;

f= [];
objPartials=[];
gbar=[];
gZ =[];
gbarZ=[];
Aact1=[];
uds = [];
lds = [];
gbarZm = {};
%% Algorithm main loop
for k = 1:opt.max_iter
    
    %%% Meanwhile condensing, study to remove this
    [xs,vs,s2,xd,vd,sd,ax,Ax,av,Av,as,As]  = condensing_R(x,u,v,s,sss,'simVars',simVars,'computeCorrection',true,'computeNullSpace',opt.condense);
    
    if imMaster
    [f,objPartials] = obj(s,u,'gradients',true);
    
    gbar.Ju = cellfun(@plus,objPartials.Ju,mudu,'UniformOutput',false);
    gbar.Js = objPartials.Js+cell2mat(muds);
    else
    gbar.Js = nan;
    gbar.Ju = zeroJacsS(1,uDims);
    objPartials.Js = nan;
    objPartials.Ju = zeroJacsS(1,uDims);
    end
    gbar.Jx = mudx;
    gbar.Jv = mudv;
    
    
    if opt.condense
        if imMaster
        gZ = mat2cell(objPartials.Js*cell2mat(As) + cell2mat(objPartials.Ju),1,uDims);
        else
        gZ = [];
        end
        
        gbarZ = cellfun(@(gx,gv)gx+gv,...
            calcgbarZ(gbar.Jx,Ax,ss,uDims,jobSchedule),...
            calcgbarZ(gbar.Jv,Av,ss,uDims,jobSchedule),...
            'UniformOutput',false);
        
        if imMaster
        gbarZ = cellfun(@(gxv,gz,gu)gxv+gz+gu,...
            gbarZ,...
            mat2cell(gbar.Js*cell2mat(As),1,uDims),...
            gbar.Ju,...
            'UniformOutput',false);
        end
        
    else
        %spmd
        objJx = zeroJacs(1,xDims);
        objJv = zeroJacs(1,vDims);
        %end
        objPartials.Jx = objJx;
        objPartials.Jv = objJv;
                
        [ sensitivities ] = generateSimulationSentivity(u,x,v,sss,simVars,[objPartials;gbar],xDims,vDims,uDims,lowActive,upActive );
        
        
        gZ = sensitivities{1};
        gbarZ = sensitivities{2};
        Aact1 = sensitivities{3};
        
        
        lowActiveSOC = lowActive;
        upActiveSOC = upActive;
        % if SOC is executed, start it with Aact1. TODO: use the jacobians
        % from the last QP, the later provides more information
        
        predictor = @(du) linearPredictor(du,x,u,v,s,sss,simVars);
        constraintBuilder = @(activeSet) generateSimulationSentivity(u,x,v,sss,simVars,[],xDims,vDims,uDims,activeSet);
        
    end
    lagFunc = @(J)simulateSystemZ_R(u,x,v,sss,J,simVars);
      
    if opt.computeCrossTerm
        
        % TODO: after finished check all input and outputs, in particular
        % uSliced!
        
        
        % Honor hard bounds in every step. Cut step if necessary
        [w,stepY] = computeCrossTerm(x,u,v,s,ax,av,as,gbarZ,sss,obj,mudx,mudu,mudv,muds,      lbx,    lbv,    lbs,    ubx,    ubv,    ubs,'xs',xs,'vs',vs,'s2',s2);
        zeta = 1;%computeZeta( gZ,M,w );
    else
        
        zeta = 1;
        stepY = 0;
        w = mat2cell(zeros(1,sum(uDims)),1,uDims);
    end
    
    
    % plot initial iterate
    if ~isempty(opt.plotFunc) && k == 1 && opt.plot
        opt.plotFunc(x,u,v,xd);
    end
    
    
    %% Update hessian approximation
    
    if k>1  % Do not perform updates if the watchdog is active!
        
        
        if imMaster
        y = cellfun(@(gbarZi,gbarZmi,wbari)gbarZi-gbarZmi-wbari,gbarZ,gbarZm,wbar,'UniformOutput',false);
        y = cell2mat(y);
        su = uV-uBV;

        % Perform the BFGS update and save information for restart
        if hInit
            M = [];
            [M,S,Y, skipping,sTy] = dampedBFGSLimRestart(M,y,su,nru,S,Y,'scale',opt.BFGSRestartscale,'it',k,'m',opt.BFGSRestartmemory,'allowDamp',opt.allowDamp);
            hInit = skipping;
        else
            [ M,S,Y,~,sTy ] = dampedBFGSLimRestart(M,y,su,nru,S,Y,'scale',opt.BFGSRestartscale,'it',k,'m',opt.BFGSRestartmemory,'allowDamp',opt.allowDamp);
        end
        end
    end
    
    
    %% Compute search direction  && lagrange multipliers
    
    % Compute bounds for the linearized problem
    udu =  cellfun(@(w,e)(w-e),ubu,u,'UniformOutput',false);
    ldu =  cellfun(@(w,e)(w-e),lbu,u,'UniformOutput',false);
    
    %spmd
    udx =  minusR(ubx,x);
    ldx =  minusR(lbx,x);
    udv =  minusR(ubv,v);
    ldv =  minusR(lbv,v);
    %end
    if imMaster
    uds =  ubs-s;
    lds =  lbs-s;
    end
    
    % Solve the QP to obtain the step on the nullspace.
    [ du,dx,dv,ds,xi,lowActive,upActive,muH,violation,qpVAl,dxN,dvN,dsN,slack,QPIT] = qpStep_R(M,gZ,w,...
        ldu,udu,...
        Aact1,predictor,constraintBuilder,...
        ax,Ax,ldx,udx,...
        av,Av,ldv,udv,...
        as,As,lds,uds,...
        sss,...
        'lowActive',lowActive,'upActive',upActive,...
        'qpDebug',opt.qpDebug,'it',k,...
        'nCons',nCons,...
        'feasTol',opt.qpFeasTol,'condense',opt.condense,'lagFunc',lagFunc,'testQP',opt.testQP);
    

    
    if imMaster && ((violation.x > masterTol) || (violation.v > masterTol) || (violation.s > masterTol)) 
        warning('QP solver too inaccurate, check the scaling and tolerance settings');
    end
    
    
    % Honor hard bounds in every step. Cut step if necessary, use the QP
    % tolerance setting to do so
    [maxStepu,du] = maximumStepLength(u,du,lbu,ubu,'tol',opt.qpFeasTol);
    
    violationx = max(violation.x,opt.qpFeasTol);
    violationv = max(violation.v,opt.qpFeasTol);
    %spmd
    [maxStepx,dx] = checkMaxStepX(x,dx,lbx,ubx,violationx);
    maxStepx = min([cell2mat(maxStepx);inf]);
    
    [maxStepv,dv] = checkMaxStepV(v,dv,lbv,ubv,violationv);
	maxStepxv = min([cell2mat(maxStepv);maxStepx]);
    
    
    maxStepxv = gopMPI('N',maxStepxv,jobSchedule);
    %end
    
    if imMaster
	violations = max(violation.s,opt.qpFeasTol);
    [maxSteps,ds] = maximumStepLength({s},{ds},{lbs},{ubs},'tol',violations);
    ds = cell2mat(ds);
    maxSteps = min(maxSteps);
    maxStep = min([maxSteps;maxStepu;maxStepxv]);
    else
    maxStep = inf;    
    end
    maxStep = NMPI_Bcast(maxStep,1,jobSchedule.Master_rank,jobSchedule.my_rank);    
    
    
    if maxStep <= 0
        printLogLine(k,...
            {'vx','vv','vs','qpFeasTol'},...
            {violation.x,violation.v,violation.s,opt.qpFeasTol},...
            'logName',['violations',num2str(jobSchedule.my_rank+1),'.txt'],'header',true);
    end
    
    %% Convergence test
    % I choose the infinity norm, because this is easier to relate to the
    % physical variables
    normdu = norm(cellfun(@(z)norm(z,'inf'),du),'inf');
    %spmd
    normax = normInf(ax);
    normav = normInf(av);
    normax = gopMPI('M',normax,jobSchedule);
    normav = gopMPI('M',normav,jobSchedule);
    %end
    normas = norm(as,'inf');
    
	converged = (normdu < opt.tolU && normax < opt.tolX && normav < opt.tolV && normas < opt.tolS  && normdu < opt.tol && normax < opt.tol && normav < opt.tol &&  normas < opt.tol && relax);
    converged = NMPI_Bcast(converged+0,1,jobSchedule.Master_rank,jobSchedule.my_rank)~=0;    
    if converged
        break;
    end
    
    %% Preparing for line-search
    
    
    
    % gbar = g+nu
    if imMaster
    gbar.Js = objPartials.Js+cell2mat(muH.ds);
    gbar.Ju = cellfun(@plus,objPartials.Ju,muH.du,'UniformOutput',false);
    else
    gbar.Js = nan;
    gbar.Ju = muH.du;
    end
    gbar.Jx = muH.dx;
    gbar.Jv = muH.dv;
    
    if relax || k == 1
        
        if  k > opt.multiplierFree
            gbarLambda.Jx = gbar.Jx;
            gbarLambda.Ju = gbar.Ju;
            gbarLambda.Jv = gbar.Jv;
            [~,~,~,lambdaX,lambdaV]= simulateSystemZ(u,x,v,sss,[],'simVars',simVars,'JacTar',gbarLambda);
            
            %{

            % first order optimality condition!
            [~,~,Jac,~,~,~] = simulateSystem(x,u,ss,'gradients',true,'xLeftSeed',lambdaX,'vLeftSeed',lambdaV,'guessX',xs,'guessV',vs,'simVars',simVars,'withAlgs',withAlgs);

            optCond.x =  cellfun(@(gbari,lambdaCi,lambdai)(gbari+(lambdaCi-lambdai)),gbarLambda.Jx,Jac.Jx,lambdaX,'UniformOutput',false);
            optCond.u =  cellfun(@(gbari,lambdaCi)(gbari+lambdaCi),gbarLambda.Ju,Jac.Ju,'UniformOutput',false);
            if withAlgs
                optCond.v = cellfun(@(gbari,lambdai)(gbari-lambdai),gbarLambda.Jv,lambdaV,'UniformOutput',false);
            end

            %}
            normInfLambda = max(cellfun(@(xv)max(abs(xv)),[lambdaX,lambdaV]));
            normInfLambda = gopMPI('M',normInfLambda,jobSchedule);
            
        else
            normInfLambda = -inf;
            
        end
        
        
        if xi ~= 1
            % multiplier free approximations
            [gbarR,errorSum,crossProduct] = multiplierFreeApproxs(gbar,ax,av,as,xd,vd,sd,w,du,xi,jobSchedule);
            % calculate equality constraints penalty
            if imMaster
            [rho,errorSumB,dualApproxB] = equalityConsPenalty(gbarR,errorSum,crossProduct,rho,rhoHat,errorSumB,dualApproxB,normInfLambda);
            else
            rho = nan;
            end
        else
            warning('xi == 1. The problem may be infeasible to solve');
        end
    end
    
    %% Merit function definition
    merit = @(f,dE,dS,varargin) l1merit(f,dE,dS,rho,jobSchedule,varargin{:});
    % line function
    phi = @(l,varargin) lineFunctionWrapper(l,...
        x,...
        v,...
        u,...
        s,...
        dx,...
        dv,...
        du,...
        ds,...
        simFunc,obj,merit,...
        jobSchedule,...
        'gradients',true,'plotFunc',opt.plotFunc,'plot',opt.plot,...
        'debug',debug,...
        'xd0',xd,...
        'vd0',vd,...
        'sd0',sd,...
        'xs0',xs,...
        'vs0',vs,...
        's20',s2,...
        'xi',xi,...
        varargin{:});
    
    
    % do not perform a watch-dog step on the very first iteration!
    if k<=1
        skipWatchDog = true;
    else
        skipWatchDog = false;
    end
    
    % Line-search
    [l,~,~,~,xfd,vars,simVars,relax,returnVars,wentBack,debugInfo] = watchdogLineSearch(phi,relax,...
        'tau',opt.tauL,'eta',opt.eta,'kmax',opt.lkMax,'debugPlot',opt.debugLS,'debug',debug&&imMaster,...
        'simVars',simVars,'curvLS',opt.curvLS,'returnVars',returnVars,'skipWatchDog',skipWatchDog,'maxStep',maxStep,'k',k);
    
    
    if relax == false  %% Watchdog step activated, should we perform SOC?
        if imMaster
        SOCCond = (debugInfo{2}.eqNorm1 > debugInfo{1}.eqNorm1);
		else
		SOCCond = 0;
        end
        SOCCond = NMPI_Bcast(SOCCond+0,1,jobSchedule.Master_rank,jobSchedule.my_rank)~=0;
        if SOCCond
        
        % build the new problem!
        varsxs = vars.xs;
        varsvs = vars.vs;
        varsx = vars.x;
        varsv = vars.v;
        %spmd
        xdSOC =  zdSOC(varsxs,varsx,xd,xi);
        vdSOC =  zdSOC(varsvs,varsv,vd,xi);
        %end
        sdSOC = zdSOCS({vars.s2},{vars.s},{sd},xi);
        sdSOC = cell2mat(sdSOC);
        
        
        [~,~,~,~,~,~,axSOC,~,avSOC,~,asSOC,~] = condensing_R(x,u,v,s,sss,...
            'simVars',simVars,...
            'computeCorrection',true,...
            'computeNullSpace',false,...
            'xd',xdSOC,'vd',vdSOC,'sd',sdSOC);
        
        
        if opt.computeCrossTerm
            [wSOC,stepYSOC] = computeCrossTerm(x,u,v,s,...
                axSOC,avSOC,asSOC,...
                gbarZ,sss,obj,...
                mudx,mudu,mudv,muds,...
                lbx,lbv,lbs,ubx,ubv,ubs,...
                'xs',xs,'vs',vs,'s2',s2);
        else
            stepYSOC = 0;
            wSOC = mat2cell(zeros(1,sum(uDims)),1,uDims);
        end
        
        
        if opt.condense
            upActiveSOC = upActive;
            lowActiveSOC = lowActive;
            Aact1 = [];
        end
        
        % Solve the QP to obtain the step on the nullspace.
        [ duSOC,dxSOC,dvSOC,dsSOC,xiSOC,lowActiveSOC,upActiveSOC,muHSOC,violationSOC,qpVAlSOC,dxNSOC,dvNSOC,dsNSOC,slack,QPITSOC] = qpStep_R(M,gZ,wSOC,...
            ldu,udu,...
            Aact1,predictor,constraintBuilder,...
            axSOC,Ax,ldx,udx,...
            avSOC,Av,ldv,udv,...
            asSOC,As,lds,uds,...
            sss,...
            'lowActive',lowActiveSOC,'upActive',upActiveSOC,...
            'qpDebug',opt.qpDebug,'it',k,...
            'nCons',nCons,...
            'feasTol',opt.qpFeasTol,'condense',opt.condense,'lagFunc',lagFunc,'testQP',opt.testQP);
        
        QPIT = QPIT+ QPITSOC;
        
        
        
        if imMaster && ((violationSOC.x > masterTol) || (violationSOC.v > masterTol) || ((violationSOC.s > masterTol)))
            warning('QP solver too inaccurate, check the scaling and tolerance settings');
        end
        
        [maxStepu,duSOC] = maximumStepLength(u,duSOC,lbu,ubu,'tol',opt.qpFeasTol); 

        % Honor hard bounds in every step. Cut step if necessary, use the QP
        % tolerance setting to do so
        violationSOCx = max(violationSOC.x,opt.qpFeasTol);
        violationSOCv = max(violationSOC.v,opt.qpFeasTol);
        %spmd
        [maxStepx,dxSOC] = checkMaxStepX(x,dxSOC,lbx,ubx,violationSOCx);
        maxStepx = min([cell2mat(maxStepx);inf]);
        
        [maxStepv,dvSOC] = checkMaxStepV(v,dvSOC,lbv,ubv,violationSOCv);
        maxStepxv = min([cell2mat(maxStepv);maxStepx]);
        
        

        maxStepxv = gopMPI('N',maxStepxv,jobSchedule);
        %end
        if imMaster
        violationSOCs = max(violationSOC.s,opt.qpFeasTol);
        [maxSteps,dsSOC] = maximumStepLength({s},{dsSOC},{lbs},{ubs},'tol',violationSOCs);
        dsSOC = cell2mat(dsSOC);
		maxSteps =  min(maxSteps);
        maxStep = min([maxSteps;maxStepu;maxStepxv]); 
        else
        maxStep = inf;
        dsSOC = nan;
        end
        maxStep = NMPI_Bcast(maxStep,1,jobSchedule.Master_rank,jobSchedule.my_rank);

        
        
        [ fSOC,dfSOC,varsSOC,simVarsSOC,debugInfoSOC ] = lineFunctionWrapper(maxStep,...
            x,...
            v,...
            u,...
            s,...
            dxSOC,...
            dvSOC,...
            duSOC,...
            dsSOC,...
            simFunc,obj,merit,...
            jobSchedule,...
            'gradients',true,'plotFunc',opt.plotFunc,'plot',opt.plot,...
            'debug',debug,...
            'xi',xi);
        
        %Try full Step
        
        xfd = [xfd;maxStep fSOC dfSOC];
        
        armijoF = @(lT,fT)  (fT - (xfd(1,2) + opt.eta*xfd(1,3)*lT));
        armijoOk = @(lT,fT) (armijoF(lT,fT) <= 0);
        
        
        debugInfoSOC.armijoVal = armijoF(maxStep,fSOC);
        debugInfo = [debugInfo;debugInfoSOC];
        
        
        if armijoOk(maxStep,fSOC)  %% accept this step!
            ax = axSOC;
            av = avSOC;
            as = asSOC;
            du = duSOC;
            dx = dxSOC;
            dv = dvSOC;
            ds = dsSOC;
            xi = xiSOC;
            lowActive = lowActiveSOC;
            upActive = upActiveSOC;
            muH = muHSOC;
            violation = violationSOC;
            qpVAl = qpVAlSOC;
            dxN = dxNSOC;
            dvN = dvNSOC;
            dsN = dsNSOC;
            w = wSOC;
            l=maxStep;
            vars = varsSOC;
            simVars = simVarsSOC;
            relax = true;
            returnVars = [];
            wentBack = false;
            
            % gbar = g+nu
            if imMaster
            gbar.Ju = cellfun(@plus,objPartials.Ju,muH.du,'UniformOutput',false);
            gbar.Js = objPartials.Js+cell2mat(muH.ds);
            else
            gbar.Js = nan;
            gbar.Ju = zeroJacsS(1,uDims);
            end
            gbar.Jx = muH.dx;
            gbar.Jv = muH.dv;
            
			if imMaster
            debugWatchdog( k,'C',xfd(end,1),xfd(end,2),xfd(end,3),debugInfo(end));
			end
        else
			if imMaster
            debugWatchdog( k,'X',xfd(end,1),xfd(end,2),xfd(end,3),debugInfo(end));
			end
        end
        
        
        end
    end
    
    
    % Restore previous lagrange multiplier estimate, if the Watch-dog
    % returned from a previous estimate
    if wentBack
        muds = muReturnS;
        mudx = muReturnX;
        mudv = muReturnV;
        mudu = muReturnU;
        muH = muHReturn;
        gbarZm = gbarZmReturn;
        gbarZ = gbarZReturn;
        
    else
        
        % calculate the lagrangian with the updated values of mu, this will
        % help to perform the BFGS update
        if opt.condense
            
            gbarZm = cellfun(@(gx,gv)gx+gv,...
                calcgbarZ(gbar.Jx,Ax,ss,uDims,jobSchedule),...
                calcgbarZ(gbar.Jv,Av,ss,uDims,jobSchedule),...
                'UniformOutput',false);
            
            if imMaster
                gbarZm = cellfun(@(gxv,gz,gu)gxv+gz+gu,...
                    gbarZm,...
                    mat2cell(gbar.Js*cell2mat(As),1,uDims),...
                    gbar.Ju,...
                    'UniformOutput',false);
            end
            
            
        else
            gbarZm = simulateSystemZ_R(u,x,v,sss,gbar,simVars);
            
        end
    end
    % Save Lagrange multipliers to restore if necessary
    if ~isempty(returnVars) && ~relax
        muReturnS = muds;
        muReturnX = mudx;
        muReturnV = mudv;
        muReturnU = mudu;
        muHReturn = muH;
        gbarZmReturn = gbarZm;
        gbarZReturn = gbarZ;
    else
        muReturnS = [];
        muReturnX = [];
        muReturnU = [];
        muReturnV = [];
        muHReturn = [];
        gbarZmReturn = [];
        gbarZReturn = [];
    end
    
    %Update dual variables estimate
    
    muHdx = muH.dx;
    muHdv = muH.dv;
    %spmd
    mudx = convexCombination(mudx,muHdx,l);
    mudv = convexCombination(mudv,muHdv,l);
    %end
    
    mudu = convexCombinationS(mudu,muH.du,l);
    if imMaster
    muds = convexCombinationS(muds,muH.ds,l);
    else
    muds = {nan};
    end
    
    %% TODO: check
    if imMaster
    gbarZm = convexCombinationS(gbarZ,gbarZm,l);
    if l == 1
        wbar = w;
    else
        wbar = cellfun(@(wi)l*wi,w,'UniformOutput',false);
    end
    else
    gbarZm = {nan};
    wbar = {nan};
    end
    
    if debug
        if imMaster
        normc = sum(dot(sd,sd));
        norma = sum(dot(as,as));
        normz = dotSumS({dsN});
        else
        normc = 0;
        norma = 0;
        normz = 0;
        end
        normc = sqrt(sum(dotSum(xd,jobSchedule)+dotSum(vd,jobSchedule)+normc));
        norma = sqrt(sum(dotSum(ax,jobSchedule)+dotSum(av,jobSchedule)+norma));
        normz = sqrt(sum(dotSum(dxN,jobSchedule)+dotSum(dvN,jobSchedule)+normz+dotSumS(du)));
        if imMaster
        printLogLine(k,...
            {'|(g+nu)Z|','|c|','|Ypy|','|Zpz|','xi','|gZ|','|w|','pz''w','stepY','l','|s|','|y|','s''y','cond(B)'},...
            {...
            sqrt(sum(cellfun(@(x)dot(x,x),gbarZm))),...
            normc,...
            norma,...
            normz,...
            xi,...
            sqrt(sum(cellfun(@(x)dot(x,x),gZ))),...
            sqrt(sum(cellfun(@(x)dot(x,x),w))),...
            sum(cellfun(@mtimes,w,du')),...
            stepY,...
            l,...
            norm(su),...
            norm(y),...
            sTy,...
            cond(M)...
            }...
            );
        end
    end
    
    % save last value of u for the BFGS update
    uBV = cell2mat(u);
    
    % return the new iterate returned after line-search.
    x = vars.x;
    xs = vars.xs;
    v =  vars.v;
    vs = vars.vs;
    if imMaster
    s =  vars.s;
    s2 = vars.s2;
    else
    s = nan;
    s2 = nan;
    end
    u = vars.u;
    
    %spmd
    x=choppBounds(lbx,x,ubx,debug);
    v=choppBounds(lbv,v,ubv,debug);
    %end
    if imMaster
    [~,s]=checkBounds(lbs ,s,ubs,'chopp',true,'verbose',debug);
    end
    
    uV = cell2mat(u);
    
    usliced = [];
    
    % Save the current iteration to a file, for debug purposes.
    if opt.saveIt
		work2Job = jobSchedule.work2Job;
        for r = 1:numel(work2Job{jobSchedule.my_rank+1})
            saveItVars(u,x{r},xs{r},v{r},vs{r},simVars{r},...
                'dir','./iterates/',...
                'it',k,...
                'r',work2Job{jobSchedule.my_rank+1}(r),...
                'keepPreviousIt',false);
        end
    end
    if ~isempty(opt.controlWriter) && imMaster
        opt.controlWriter(u,k);
    end
    
    
    % print main debug
    if  debug && imMaster
        if mod(k,10) == 1
            header = true;
        else
            header = false;
        end
        tMax = 0;
        
        violationH = violation.x;
        violationH = max(violationH,violation.v);
        
        dispFunc(k,norm(cell2mat(gbarZm)),violationH,normdu,rho,tMax,xfd,cond(M),relax,debugInfo,header,QPIT );
    end
    
    if imMaster && l == 0  %line search couldn't make a sufficient decrease
        warning('lineSearch determined 0 step length');
        break;
    end
    
    
end


% recover previous variables if you performed a watch-dog step in the last iteration
if ~converged &&  ~relax
    x = returnVars.vars0.x;
    u = returnVars.vars0.u;
    v = returnVars.vars0.v;
    s = returnVars.vars0.s;
    
    simVars = returnVars.simVars0;
    [f] = obj(s,u);
    
end

end

function me = convexCombination(m,mH,l)
f = @(mr,mHr)convexCombinationS(mr,mHr,l);
me = cellfun(f,m,mH,'UniformOutput',false);
end

function me = convexCombinationS(m,mH,l)
me = cellfun(@(x1,x2)(1-l)*x1+l*x2,m,mH,'UniformOutput',false);
end


function zm = minusR(z1,z2)
    if ~isempty(z2)
        if isnumeric(z1{1})
            zm = cellfun(@(z2)minusS(z1,z2),z2,'UniformOutput',false);
        else
            zm = cellfun(@minusS,z1,z2,'UniformOutput',false);
        end
    else
        zm = cell(0,1);
    end
end

function me = minusS(mU,mL)
me = cellfun(@minus,mU,mL,'UniformOutput',false);
end

function me = calcgbarZ(J,A,ss,uDims,jobSchedule)
%spmd

    f = @calcgbarZS;
    vr = cellfun(f,J,A,ss,'UniformOutput',false);
    me = catAndSum(vr);
    me = gopMPI('+',me,jobSchedule);
%end

me = mat2cell(me,size(me,1),uDims);
end
function v = calcgbarZS(J,A,ss)
v = cellmtimesT( J,A,'lowerTriangular',true,'ci',ss.ci,'columnVector',false);
end

function out = catAndSum(M)
if ~isempty(M)
    M = cellfun(@cell2mat,M,'UniformOutput',false);
    
    if any(cellfun(@issparse,M))
        if isrow(M)
            M = M';
        end
        rows= size(M{1},1);
        blocks = numel(M);
        out = sparse( repmat(1:rows,1,blocks),1:rows*blocks,1)*cell2mat(M);
    else
        out = sum(cat(3,M{:}),3);
    end
    
else
    out = 0;
end
end


function e = dotSum(z,jobSchedule)

%spmd
f = @dotSumS;
e = sum(cellfun(f,z));
e = gopMPI('+',e,jobSchedule);
%end

end
function e = dotSumS(z)
e = sum(cellfun(@(zi)sum(dot(zi,zi)),z));
end


function zd = zdSOC(lzs,lz,zd,xi)
f= @(lzs,lz,zd)zdSOCS(lzs,lz,zd,xi);
zd = cellfun(f,lzs,lz,zd,'UniformOutput',false);
end
function zd = zdSOCS(lzs,lz,zd,xi)
zd = cellfun(@(lzsi,lzi,zdi)lzsi-lzi+(1-xi)*zdi,lzs,lz,zd,'UniformOutput',false);
end

function xDims = getXDims0(ss)
    xDims = cellfun(@(ssr)numel(ssr.state),ss); 
end

function zb = buildRelaxedBound(zDims,sign)
zb = cellfun(@(zDimsr)buildRelaxedBoundS(zDimsr,sign),zDims,'UniformOutput',false);
end
function zb = buildRelaxedBoundS(zDims,sign)
zb = arrayfun(@(nzk)sign*inf(nzk,1),zDims,'UniformOutput',false);
end

function zDims = getZDims(z)
zDims = cellfun(@(zr)cellfun(@numel,zr),z,'UniformOutput',false);
end

function [z,ok] = choppBounds(lbz,z,ubz,debug)
if ~isempty(lbz)
    if isnumeric(lbz{1})
        [ok,z]  = cellfun(@(zr)        checkBounds(lbz,zr,ubz,'chopp',true,'verbose',debug),    z,    'UniformOutput',false);
        ok = cellfun(@all,ok);
        ok = all(ok);
    else % must be a cell
        [ok,z]  = cellfun(@(lbr,zr,ubr)checkBounds(lbr,zr,ubr,'chopp',true,'verbose',debug),lbz,z,ubz,'UniformOutput',false);
        ok = cellfun(@all,ok);
        ok = all(ok);
    end
else
    z = cell(0,1);
    ok = true;
end
end

function zAct = startEmptyActiveSet(zDims)
zAct = cellfun(@(zDimsr)arrayfun(@(zrkn)false(zrkn,1),zDimsr,'UniformOutput',false),zDims,'UniformOutput',false);
end

function mudz = initDualVariable(zDims)
    mudz = cellfun(@initDualVariableS,zDims,'UniformOutput',false);
end
function mudz = initDualVariableS(zDimsr)
    mudz = sparse(1,sum(zDimsr));
    mudz = mat2cell(mudz,1,zDimsr);
end

function [maxStepx,dx] = checkMaxStepX(x,dx,lbx,ubx,violationx)
    [maxStepx,dx] = cellfun(@(zi,dz)maximumStepLength(zi,dz,lbx,ubx,'tol',violationx),x,dx,'UniformOutput',false);
end

function [maxStepv,dv] = checkMaxStepV(v,dv,lbv,ubv,violationv)
	[maxStepv,dv] = cellfun(@(zi,dz,lb,ub)maximumStepLength(zi,dz,lb,ub,'tol',violationv),v,dv,lbv,ubv,'UniformOutput',false);
end

function normz = normInf(z)
    normz = norm([cellfun(@(z)norm(cell2mat(z),'inf'),z);0],'inf');
end

function J = zeroJacs(n,zDims)
f = @(zDimsr)zeroJacsS(n,zDimsr);
J = cellfun(f,zDims,'UniformOutput',false);
end
function J = zeroJacsS(n,zDims)
dimsSum = sum(zDims);
    J = sparse(n,dimsSum);
J = mat2cell(J,n,zDims);
    
end
