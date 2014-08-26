function [u,x,v,f,xd,M,simVars] = remso(u,ss,obj,varargin)
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
%   lbxH - State hard lower bound for each point in the prediction horizon.
%
%   ubxH - State hard  upper bound for each point in the prediction horizon.
%
%   lbvH - Algebraic state hard lower bound for each point in the prediction horizon.
%
%   ubvH - Algebraic state hard upper bound for each point in the prediction horizon.
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
             'lbxH',[],'ubxH',[],'lbvH',[],'ubvH',[],...
    'tol',1e-1,'tolU',1e-2,'tolX',1e-2,'tolV',1e-2,'max_iter',50,...
    'M',[],'x',[],'v',[],...
    'plotFunc',[],...
    'lkMax',4,'eta',0.1,'tauL',0.1,'debugLS',false,'curvLS',true,...
    'qpDebug',true,...
    'lowActive',[],'upActive',[],...
    'simVars',[],'debug',true,'plot',false,'saveIt',false);

opt = merge_options(opt, varargin{:});


jobSchedule = ss.jobSchedule;

% extract information on the prediction horizon and control intervals
totalPredictionSteps = getTotalPredictionSteps(ss);
totalControlSteps = numel(u);

% number of variables
nx = numel(ss.state);
nu = numel(u{1});
nv = ss.nv;

% true if dealing with algebraic states
withAlgs = (nv>0);

% dimension of the control space, dimension of the reduced problem
nru = numel(cat(2,u{:}));

%% Control, state and algebraic state bounds processing
uV = cell2mat(u);
if isempty(opt.lbu)
    lbu = [];
else
    lbu = cell2mat(opt.lbu);
    if ~all(uV-lbu >=0)
        error('Make a feasible first guess of the control variables')
    end
end
if isempty(opt.ubu)
    ubu = [];
else
    ubu = cell2mat(opt.ubu);
    if ~all(ubu-uV >=0)
        error('Make a feasible first guess of the control variables')
    end
end
if isempty(opt.lbx)
    opt.lbx = repmat({-inf(nx,1)},totalPredictionSteps,1);
end
if isempty(opt.ubx)
    opt.ubx = repmat({inf(nx,1)},totalPredictionSteps,1);
end
if withAlgs && isempty(opt.lbv)
    opt.lbv = repmat({-inf(nv,1)},totalPredictionSteps,1);
end
if withAlgs && isempty(opt.ubv)
    opt.ubv = repmat({inf(nv,1)},totalPredictionSteps,1);
end

checkHardConstraints = false;
if isempty(opt.lbxH)
    opt.lbxH = repmat({-inf(nx,1)},totalPredictionSteps,1);
else
    checkHardConstraints = true;
end
if isempty(opt.ubxH)
    opt.ubxH = repmat({inf(nx,1)},totalPredictionSteps,1);
else
    checkHardConstraints = true;    
end
if withAlgs && isempty(opt.lbvH)
    opt.lbvH = repmat({-inf(nv,1)},totalPredictionSteps,1);
else
    checkHardConstraints = true;    
end
if withAlgs && isempty(opt.ubvH)
    opt.ubvH = repmat({inf(nv,1)},totalPredictionSteps,1);
else
    checkHardConstraints = true;
end

% solf bounds must be bounded by hard bounds
if checkHardConstraints
    
    opt.lbx = cellfun(@(l1,l2)max(l1,l2),opt.lbx,opt.lbxH,'UniformOutput',false);
    opt.ubx = cellfun(@(l1,l2)min(l1,l2),opt.ubx,opt.ubxH,'UniformOutput',false);
	
	if withAlgs
        opt.lbv = cellfun(@(l1,l2)max(l1,l2),opt.lbv,opt.lbvH,'UniformOutput',false);
        opt.ubv = cellfun(@(l1,l2)min(l1,l2),opt.ubv,opt.ubvH,'UniformOutput',false);
	end
    
end


maxStep = 1;

udv = [];
ldv = [];
dv = [];

% Multiple shooting simulation function
simFunc = @(xk,uk,varargin) simulateSystem(xk,uk,ss,varargin{:});


%% Define empty active sets if they are not given
if isempty(opt.lowActive)
    opt.lowActive.x = cellfun(@(x)false(size(x)),opt.lbx,'UniformOutput',false);
    if withAlgs
        opt.lowActive.v = cellfun(@(x)false(size(x)),opt.lbv,'UniformOutput',false);
    end
end
if isempty(opt.upActive)
    opt.upActive.x = cellfun(@(x)false(size(x)),opt.ubx,'UniformOutput',false);
    if withAlgs
        opt.upActive.v = cellfun(@(x)false(size(x)),opt.ubv,'UniformOutput',false);
    end
end

%% initial simulation profile
if isempty(opt.simVars)
    simVars = cell(totalPredictionSteps,1);
else
    simVars = opt.simVars;
end

%% Process initial MS simulation guess, if not given, get it by forward simulation
simulateSS = false;
if ~isempty(opt.x)
    %  Initial guess for prediction given by the user
    x = opt.x;
    xs.client = opt.x;
else
    % Initial guess not provided, take from a simulation in the gradient
    % routine
    simulateSS = true;
    xs.client = cell(totalPredictionSteps,1);
end
if withAlgs
    if isempty(opt.v)
        v = repmat({zeros(nv,1)},totalPredictionSteps,1);
        vs.client = cell(totalPredictionSteps,1);
    else
        v = opt.v;
        vs.client = opt.v;
    end
else
    v = [];
    vs = [];
end

if simulateSS
	[~,~,~,simVars,xsR,vsR,uslicedR] = simulateSystemSS(u,ss,[],'guessX',xs.client,'guessV',vs.client,'simVars',simVars);
    x = xsR;
    v = vsR;
    xs.client = xsR;
    vs.client = vsR;
    usliced.client = uslicedR;
else
    [xsR,vsR,~,~,simVars,uslicedR] = simulateSystem(x,u,ss,'gradients',false,'guessX',xs.client,'guessV',vs.client,'simVars',simVars);
	xs.worker = xsR;
    vs.worker = vsR;
    xs = rmfield(xs,'client');
	vs = rmfield(vs,'client');
    usliced.worker = uslicedR;
end


%% lagrange multipliers estimate initilization 

mudx= repmat({zeros(nx,1)},totalPredictionSteps,1);
mudu = repmat({zeros(nu,1)},totalControlSteps,1);
if withAlgs
    mudv = repmat({zeros(nv,1)},totalPredictionSteps,1);
end


%% Hessian Initializaiton
if(isempty(opt.M))
    hInit = true;
    M = eye(nru);
else
    hInit = false;
    M = opt.M;
end
%% Curvature history record
S = [];
Y = [];



%% Line-search parameters
rho = 1/(totalPredictionSteps*(nx+nv));
rhoHat = rho/100;
returnVars = [];
relax = false;   % to avoid the hessian update and perform a fine line-search
errorSumB = [];
dualApproxB = [];
tau = [];


%%  This file allows you to stop the algorithm for debug during execution.
% If the file is deleted, the algorithm will stop at the predefined set
% points.
if opt.debug
    fid = fopen('deleteMe2Break.txt','w');fclose(fid);
end

% convergence flag
converged = false;


%% Algorithm main loop
for k = 1:opt.max_iter
    
    % Perform the condensing thechnique on the current iterate
    [xs.client,vs.client,xd,vd,ax,Ax,av,Av]  = condensing(x,u,v,ss,'simVars',simVars);

    % Calculate the objective function gradient
    [f,B,objPartials] = targetGrad(xs,u,vs,obj,Ax,Av,ss.ci,'usliced',usliced);
    
    % plot initial iterate
    if ~isempty(opt.plotFunc) && k == 1 && opt.plot
        opt.plotFunc(x,u,v,xd);
    end
    
    % debug cheack-point, check if the file is present
    if opt.debug
        fid = fopen('deleteMe2Break.txt','r');
        if fid == -1
            fid = fopen('deleteMe2Break.txt','w');fclose(fid);
            keyboard;
        else
            fclose(fid);
        end
    end
    
    %% Update hessian approximation
    
    if relax  % Do not perform updates if the watchdog is active!
        
        % Calculate the lagrangian
        L = cat(2,B{:});
        lagTauX = cellmtimesT( mudx,Ax,'lowerTriangular',true,'ci',ss.ci);
        L = L + cat(2,lagTauX{:});
        if withAlgs
            lagTauV = cellmtimesT( mudv,Av,'lowerTriangular',true,'ci',ss.ci);
            L = L + cat(2,lagTauV{:});
        end
        
        % Perform the BFGS update and save information for restart
        if hInit
            [  M,S,Y, skipping ] = dampedBFGSLimRestart(M,L-LB,uV-uBV,nru,S,Y,'scale',true,'it',k);
            hInit = skipping;
        else
            [  M,S,Y ] = dampedBFGSLimRestart(M,L-LB,uV-uBV,nru,S,Y,'scale',false,'it',k);
        end
        %        [ M,S,Y ] = limitedBFGS(L-LB,uV-uBV,nru,S,Y);
        %        [ M,S,Y ] = limitedDumpedBFGS(M,L-LB,uV-uBV,nru,S,Y );
        
        
    end
    
    
    %% Compute search direction  && lagrange multipliers
    
    % Compute bounds for the linearized problem
    udx =  cellfun(@(w,e,r)(w-e-r),opt.ubx,x,ax,'UniformOutput',false);
    ldx =  cellfun(@(w,e,r)(w-e-r),opt.lbx,x,ax,'UniformOutput',false);
    if withAlgs
        udv =  cellfun(@(w,e,r)(w-e-r),opt.ubv,v,av,'UniformOutput',false);
        ldv =  cellfun(@(w,e,r)(w-e-r),opt.lbv,v,av,'UniformOutput',false);
    end
    
    % Solve the QP to obtain the step on the nullspace.
    [ duN,dxN,dvN,opt.lowActive,opt.upActive,muH,s,violationH,qpVAl] = prsqpStep(M,B,...
        u,lbu,ubu,...
        Ax,ldx,udx,...
        Av,ldv,udv,...
        'lowActive',opt.lowActive,'upActive',opt.upActive,...
        'ci',ss.ci,...
        'qpDebug',opt.qpDebug,'it',k);
    
    % debug cheack-point, check if the file is present
    if opt.debug
        fid = fopen('deleteMe2Break.txt','r');
        if fid == -1
            fid = fopen('deleteMe2Break.txt','w');fclose(fid);
            keyboard;
        else
            fclose(fid);
        end
    end
    
    % define the PRSQP step by adding the range space solution and
    % nullspace solution
    du = duN;
    dx = cellfun(@(z,dz)z+dz,ax,dxN,'UniformOutput',false);
    if withAlgs
        dv = cellfun(@(z,dz)z+dz,av,dvN,'UniformOutput',false);
    end
    
	% Honor hard bounds in every step. Cut step if necessary
    if checkHardConstraints
        [maxStep,dx] = maximumStepLength(x,dx,opt.lbxH,opt.ubxH);
        if withAlgs
            [maxStepv,dv] =maximumStepLength(v,dv,opt.lbvH,opt.ubvH);
            maxStep = min(maxStep,maxStepv);
        end
    end
    
    %%% test  firstOrderOpt < tol !
    % [firstOrderOpt] = testFirstOrderOpt(M,objPartials,duN,dxN,dvN,muH,withAlgs)
    
    
    %% Convergence test
    % I choose the infinity norm, because this is easier to relate to the
    % physical variables
    normdu = norm(cellfun(@(z)norm(z,'inf'),duN),'inf');
    normax = norm(cellfun(@(z)norm(z,'inf'),ax),'inf');
    normav = norm(cellfun(@(z)norm(z,'inf'),av),'inf');
    
    if normdu < opt.tolU && normax < opt.tolX && normav < opt.tolV && normdu < opt.tol && normax < opt.tol && normav < opt.tol && violationH(end) < opt.tol && relax
        converged = true;
        break;
    end
    
    %% Preparing for line-search
    
    if relax || k == 1
        % multiplier free approximations
        [firstOptDualApprox,errorSum] = multiplierFreeApproxs(objPartials,ax,av,xd,vd,muH,withAlgs);
        % calculate equality constraints penalty
        [rho,errorSumB,dualApproxB] = equalityConsPenalty(firstOptDualApprox,errorSum,rho,rhoHat,errorSumB,dualApproxB);        
        % Calculate penalties for the bounds violations
        [tau] = boundViolationWeights(muH,tau,withAlgs);
        
        % Adapt bounds according to the least-infeasibility
        [lbs,ubs] = leastInfeasibleBounds(s,opt,withAlgs);
    end
    
    %% Merit function definition
    lbsD = cellfun(@(z) distributeVariables( z,jobSchedule),lbs,'UniformOutput',false);
    ubsD = cellfun(@(z) distributeVariables( z,jobSchedule),ubs,'UniformOutput',false);
    tauD = cellfun(@(z) distributeVariables( z,jobSchedule),tau,'UniformOutput',false);
    merit = @(f,dE,bE,varargin) l1merit(f,dE,bE,ubsD,lbsD,rho,tauD,varargin{:});

	xW  = distributeVariables(x ,jobSchedule);
	vW  = distributeVariables(v ,jobSchedule);
	dxW = distributeVariables(dx,jobSchedule);
	dvW = distributeVariables(dv,jobSchedule);
	xdW = distributeVariables(xd,jobSchedule);
	vdW = distributeVariables(vd,jobSchedule);
	xsW = distributeVariables(xs,jobSchedule);
	vsW = distributeVariables(vs,jobSchedule);
    
    % line function
    phi = @(l,varargin) lineFunctionWrapper(l,...
        xW,...
        vW,...
        u,...
        dxW,...
        dvW,...
        du,...
        simFunc,obj,merit,jobSchedule,'gradients',true,'plotFunc',opt.plotFunc,'plot',opt.plot,...
        'debug',opt.debug,...
        'xd0',xdW,...
        'vd0',vdW,...
        'xs0',xsW,...
        'vs0',vsW,...
        varargin{:});
   
    
    % do not perform a watch-dog step on the very first iteration! 
    if k<=1
        skipWatchDog = true;
    else
        skipWatchDog = false;
    end
    
    % Line-search 
    [l,~,~,~,xfd,vars,simVars,relax,returnVars,wentBack,debugInfo] = watchdogLineSearch(phi,relax,...
        'tau',opt.tauL,'eta',opt.eta,'kmax',opt.lkMax,'debug',opt.debugLS,...
        'simVars',simVars,'curvLS',opt.curvLS,'returnVars',returnVars,'skipWatchDog',skipWatchDog,'maxStep',maxStep);

    % debug cheack-point, check if the file is present
    if opt.debug
        fid = fopen('deleteMe2Break.txt','r');
        if fid == -1
            fid = fopen('deleteMe2Break.txt','w');fclose(fid);
            keyboard;
        else
            fclose(fid);
        end
    end
    
    % Restore previous lagrange multiplier estimate, if the Watch-dog
    % returned from a previous estimate
    if wentBack 
        mudx = muReturnX;
        if withAlgs
            mudv = muReturnV;
        end
        mudu = muReturnU;
        muH = muHReturn;
    end
    % Save Lagrange multipliers to restore if necessary
    if ~isempty(returnVars)
        muReturnX = mudx;
        if withAlgs
            muReturnV = mudv;
        end
        muReturnU = mudu;
        muHReturn = muH;
    else
        muReturnX = [];
        muReturnU = [];
        muReturnV = [];
        muHReturn = [];
    end
    
    %Update dual variables estimate
    mudx = cellfun(@(x1,x2,x3)(1-l)*x1+l*(x2-x3),mudx,muH.ub.x,muH.lb.x,'UniformOutput',false);
    mudu = cellfun(@(x1,x2,x3)(1-l)*x1+l*(x2-x3),mudu,muH.ub.u,muH.lb.u,'UniformOutput',false);  
    if withAlgs
        mudv = cellfun(@(x1,x2,x3)(1-l)*x1+l*(x2-x3),mudv,muH.ub.v,muH.lb.v,'UniformOutput',false);
    end
    
    
    % calculate the lagrangian with the updated values of mu, this will
    % help to perform the BFGS update
    LB = cat(2,B{:});
    lagTauX = cellmtimesT( mudx,Ax,'lowerTriangular',true,'ci',ss.ci);
    LB = LB + cat(2,lagTauX{:});
    if withAlgs
        lagTauV = cellmtimesT( mudv,Av,'lowerTriangular',true,'ci',ss.ci);
        LB = LB + cat(2,lagTauV{:});
    end
    
    % save last value of u for the BFGS update
    uBV = cell2mat(u);
    
    % return the new iterate returned after line-search.
    x = bringVariables(vars.x,jobSchedule);  
    xs.worker = vars.xs;
    xs = rmfield(xs,'client');
    if withAlgs
        v =  bringVariables(vars.v,jobSchedule);
        vs.worker = vars.vs;
        vs = rmfield(vs,'client');
    end
    u = vars.u;
    
    uV = cell2mat(u);
    if ~all(uV-lbu >=-eps) || ~all(ubu-uV >=-eps)
        warning('Control values out of feasible set')
    end
    usliced = vars.usliced;
    
    % Save the current iteration to a file, for debug purposes.
    if opt.saveIt
        save itVars x u v xd rho M tau;
    end
    
    % print main debug
    if  opt.debug
        if mod(k,10) == 1
            header = true;
        else
            header = false;
        end
        tMax = norm(cellfun(@(x)norm(x,'inf'),cat(2,tau{:})),'inf');
        
        L = cat(2,B{:});
        lagTauX = cellmtimesT( mudx,Ax,'lowerTriangular',true,'ci',ss.ci);
        L = L + cat(2,lagTauX{:});
        if withAlgs
            lagTauV = cellmtimesT( mudv,Av,'lowerTriangular',true,'ci',ss.ci);
            L = L + cat(2,lagTauV{:});
        end
        L = L + (cat(1,mudu{:}))';
        dispFunc(k,norm(L),violationH,normdu,rho,tMax,xfd,cond(M),relax,debugInfo,header );
    end
    
    if l == 0  %line search couldn't make a sufficient decrease
        warning('lineSearch determined 0 step length');
        break;
    end
    
    
end


% recover previous variables if you performed a watch-dog step in the last iteration
if ~converged &&  ~relax
        x = bringVariables(returnVars.vars.x,jobSchedule);
        u = returnVars.vars.u;
        if withAlgs
            v = bringVariables(returnVars.vars.v,jobSchedule);
        end
        simVars = returnVars.simVars;
        [xs,vs,~,~,simVars] = simulateSystem(x,u,ss,'guessV',v,'simVars',simVars);
        f = obj(xs,u,v,'gradients',false);
        xsF = bringVariables(xs,jobSchedule);
        xd = cellfun(@(x1,x2)x1-x2,xsF,x,'UniformOutput',false);
end


end


