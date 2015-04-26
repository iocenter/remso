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
    'BFGSRestartscale', true,'BFGSRestartmemory',6,...
    'lkMax',4,'eta',0.1,'tauL',0.1,'debugLS',false,'curvLS',true,...
    'qpDebug',true,...
    'lowActive',[],'upActive',[],...
    'simVars',[],'debug',true,'plot',false,'saveIt',false,...
    'condensingParallel',false,...
    'controlWriter',[],...
    'multiplierFree',inf,...
    'allowDamp',true,...
    'qpFeasTol',1e-6,...
    'condense',false,...
    'computeCrossTerm',true);

opt = merge_options(opt, varargin{:});


masterTol = min([opt.tol,opt.tolU,opt.tolX,opt.tolV]);

%The qpFeasTol must be tighter than tol, tolX, tolV, and tolU'
if opt.qpFeasTol > masterTol
    opt.qpFeasTol = masterTol;
end

jobSchedule = ss.jobSchedule;

% extract information on the prediction horizon and control intervals
totalPredictionSteps = getTotalPredictionSteps(ss);


% number of variables
nx = numel(ss.state);
uDims = cellfun(@(uu)size(uu,1),u);

% dimension of the control space, dimension of the reduced problem
nru = sum(uDims);

%% Control and state bounds processing
if isempty(opt.lbu)
    opt.lbu = cellfun(@(z)-inf(size(z)),u,'UniformOutput',false);
end
if isempty(opt.ubu)
    opt.ubu = cellfun(@(z)inf(size(z)),u,'UniformOutput',false);
end

[~,u]  = checkBounds( opt.lbu,u,opt.ubu,'chopp',true,'verbose',opt.debug);
uV = cell2mat(u);

if isempty(opt.lbx)
    opt.lbx = repmat({-inf(nx,1)},totalPredictionSteps,1);
end
if isempty(opt.ubx)
    opt.ubx = repmat({inf(nx,1)},totalPredictionSteps,1);
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
    x = cell(totalPredictionSteps,1);
    xs.client = cell(totalPredictionSteps,1);
end
    if isempty(opt.v)
        vs.client = cell(totalPredictionSteps,1);
    else
        vs.client = opt.v;
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
    v = bringVariables(vsR,jobSchedule);
    usliced.worker = uslicedR;
end

xDims = cellfun(@numel,x);
vDims = cellfun(@numel,v);
withAlgs = sum(vDims)>0;



%% algebraic state bounds processing
if withAlgs && isempty(opt.lbv)
    opt.lbv = arrayfun(@(d)-inf(d,1),vDims,'UniformOutput',false);
end
if withAlgs && isempty(opt.ubv)
    opt.ubv = arrayfun(@(d)inf(d,1),vDims,'UniformOutput',false);
end

[~,x]  = checkBounds( opt.lbx,x,opt.ubx,'chopp',true,'verbose',opt.debug);
if withAlgs
    [~,v]  = checkBounds( opt.lbv,v,opt.ubv,'chopp',true,'verbose',opt.debug);
end

%% hard constraints
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
    opt.lbvH = arrayfun(@(d)-inf(d,1),vDims,'UniformOutput',false);
else
    checkHardConstraints = true;    
end
if withAlgs && isempty(opt.ubvH)
    opt.ubvH = arrayfun(@(d)inf(d,1),vDims,'UniformOutput',false);
else
    checkHardConstraints = true;
end

% solve bounds must be bounded by hard bounds
if checkHardConstraints
    
    opt.lbx = cellfun(@(l1,l2)max(l1,l2),opt.lbx,opt.lbxH,'UniformOutput',false);
    opt.ubx = cellfun(@(l1,l2)min(l1,l2),opt.ubx,opt.ubxH,'UniformOutput',false);
	
	if withAlgs
        opt.lbv = cellfun(@(l1,l2)max(l1,l2),opt.lbv,opt.lbvH,'UniformOutput',false);
        opt.ubv = cellfun(@(l1,l2)min(l1,l2),opt.ubv,opt.ubvH,'UniformOutput',false);
	end
    
end


udv = [];
ldv = [];
dv = [];
Aact1= [];
predictor = [];
constraintBuilder = [];

% Multiple shooting simulation function
simFunc = @(xk,uk,varargin) simulateSystem(xk,uk,ss,'withAlgs',withAlgs,varargin{:});


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



%% lagrange multipliers estimate initilization
mudx= repmat({zeros(nx,1)},totalPredictionSteps,1);
mudu = cellfun(@(z)zeros(size(z)),u,'UniformOutput',false);
if withAlgs
    mudv = cellfun(@(z)zeros(size(z)),v,'UniformOutput',false);
else
    mudv = [];
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
if opt.debug
    fid = fopen('logBFGS.txt','w');
    fclose(fid); 
end


%% Curvature history record
y = zeros(1,sum(uDims));
s = zeros(sum(uDims),1);
sTy = 0;

S = [];
Y = [];



%% Line-search parameters
rho = 1/(totalPredictionSteps*nx+sum(vDims));
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
    
    %%% Meanwhile condensing, study to remove this
    % Perform the condensing thechnique on the current iterate

	[xs.client,vs.client,xd,vd,ax,Ax,av,Av,simVars]  = condensing(x,u,v,ss,'simVars',simVars,'computeCorrection',true,'withAlgs',withAlgs,'computeNullSpace',opt.condense);


    [f,objPartials] = obj(x,u,v,'gradients',true);
    
    gbar.Jx =  cellfun(@(Jz,m)(Jz+m'),objPartials.Jx,mudx','UniformOutput',false);
    gbar.Ju =  cellfun(@(Jz,m)(Jz+m'),objPartials.Ju,mudu','UniformOutput',false);
    if withAlgs
        gbar.Jv = cellfun(@(Jz,m)(Jz+m'),objPartials.Jv,mudv','UniformOutput',false);
    end    
    
    if opt.condense
    
        if withAlgs
            gZ = vectorTimesZ(objPartials.Jx,objPartials.Ju,objPartials.Jv,Ax,Av,ss.ci );
        else
            gZ = vectorTimesZ(objPartials.Jx,objPartials.Ju,[],Ax,[],ss.ci );
        end  
    
        if withAlgs
            gbarZ = vectorTimesZ(gbar.Jx,gbar.Ju,gbar.Jv,Ax,Av,ss.ci );
        else
            gbarZ = vectorTimesZ(gbar.Jx,gbar.Ju,[],Ax,[],ss.ci );
        end
        
    else
        
        
        [ sensitivities ] = generateSimulationSentivity(u,x,v,ss,simVars,[objPartials;gbar],xDims,vDims,uDims,opt.lowActive,opt.upActive );
        
        
        gZ = sensitivities{1};
        gbarZ = sensitivities{2};
        Aact1 = sensitivities{3};
       
        
        lowActiveSOC = opt.lowActive;
        upActiveSOC = opt.upActive;
        % if SOC is executed, start it with Aact1. TODO: use the jacobians
        % from the last QP, the later provides more information
        
    	predictor = @(du) linearPredictor(du,x,u,v,ss,simVars,withAlgs);
        constraintBuilder = @(activeSet) linearizeActiveConstraint(activeSet,u,x,v,ss,simVars,withAlgs );    
    end
    
    % TODO: after finished check all input and outputs, in particular
    % uSliced!

    
    % Honor hard bounds in every step. Cut step if necessary
    if opt.computeCrossTerm
        [w,stepY] = computeCrossTerm(x,u,v,ax,av,gbarZ,ss,obj,mudx,mudu,mudv,opt.lbxH,opt.lbvH,opt.ubxH,opt.ubvH,withAlgs,'xs',xs.client,'vs',vs.client);
        zeta = 1;%computeZeta( gZ,M,w );
    else
        w = cellfun(@(xx)zeros([size(xx,2),size(xx,1)]),u','UniformOutput',false);
        stepY = 0;       
    end

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
    
    if k>1  % Do not perform updates if the watchdog is active!
        
        
        
        y = cellfun(@(gbarZi,gbarZmi,wbari)gbarZi-gbarZmi-wbari,gbarZ,gbarZm,wbar,'UniformOutput',false);
        y = cell2mat(y);
        s = uV-uBV;
        
        % Perform the BFGS update and save information for restart
        if hInit
            M = [];
            [M,S,Y, skipping,sTy] = dampedBFGSLimRestart(M,y,s,nru,S,Y,'scale',opt.BFGSRestartscale,'it',k,'m',opt.BFGSRestartmemory,'allowDamp',opt.allowDamp);
            hInit = skipping;
        else
            [ M,S,Y,~,sTy ] = dampedBFGSLimRestart(M,y,s,nru,S,Y,'scale',opt.BFGSRestartscale,'it',k,'m',opt.BFGSRestartmemory,'allowDamp',opt.allowDamp);
        end
        
    end
    
    
    %% Compute search direction  && lagrange multipliers
    
    % Compute bounds for the linearized problem
    udu =  cellfun(@(w,e)(w-e),opt.ubu,u,'UniformOutput',false);
    ldu =  cellfun(@(w,e)(w-e),opt.lbu,u,'UniformOutput',false);
    
    udx =  cellfun(@(w,e)(w-e),opt.ubx,x,'UniformOutput',false);
    ldx =  cellfun(@(w,e)(w-e),opt.lbx,x,'UniformOutput',false);
    if withAlgs
        udv =  cellfun(@(w,e)(w-e),opt.ubv,v,'UniformOutput',false);
        ldv =  cellfun(@(w,e)(w-e),opt.lbv,v,'UniformOutput',false);
    end
    
        
    % Solve the QP to obtain the step on the nullspace.
    [ du,dx,dv,xi,opt.lowActive,opt.upActive,muH,violation,qpVAl,dxN,dvN,~,QPIT] = qpStep(M,gZ,w,...
        ldu,udu,...
        Aact1,predictor,constraintBuilder,...
        ax,Ax,ldx,udx,...
        av,Av,ldv,udv,...
        'lowActive',opt.lowActive,'upActive',opt.upActive,...
        'ci',ss.ci,...
        'qpDebug',opt.qpDebug,'it',k,'withAlgs',withAlgs,...
        'feasTol',opt.qpFeasTol,'condense',opt.condense);
    
    % debug check-point, check if the file is present
    if opt.debug
        fid = fopen('deleteMe2Break.txt','r');
        if fid == -1
            fid = fopen('deleteMe2Break.txt','w');fclose(fid);
            keyboard;
        else
            fclose(fid);
        end
    end
    
    if violation.x > masterTol || (withAlgs && (violation.v > masterTol))
        warning('QP solver too inaccurate, check the scaling and tolerance settings');
    end

    
    % Honor hard bounds in every step. Cut step if necessary, use the QP
    % tolerance setting to do so
    [maxStep,du] = maximumStepLength(u,du,opt.lbu,opt.ubu,'tol',opt.qpFeasTol);
    
    [maxStepx,dx] = maximumStepLength(x,dx,opt.lbx,opt.ubx,'tol',violation.x);
    maxStep = min(maxStep,maxStepx);
    if withAlgs
        [maxStepv,dv] =maximumStepLength(v,dv,opt.lbv,opt.ubv,'tol',violation.v);
        maxStep = min(maxStep,maxStepv);
    end
    
    
    
    
    %% Convergence test
    % I choose the infinity norm, because this is easier to relate to the
    % physical variables
    normdu = norm(cellfun(@(z)norm(z,'inf'),du),'inf');
    normax = norm(cellfun(@(z)norm(z,'inf'),ax),'inf');
    normav = norm(cellfun(@(z)norm(z,'inf'),av),'inf');
    
    if normdu < opt.tolU && normax < opt.tolX && normav < opt.tolV && normdu < opt.tol && normax < opt.tol && normav < opt.tol && relax
        converged = true;
        break;
    end
    
    %% Preparing for line-search
    
    % gbar = g+nu
    gbar.Jx =  cellfun(@(Jz,mub,mul)(Jz+(mub-mul)'),objPartials.Jx,muH.ub.x',muH.lb.x','UniformOutput',false);
    gbar.Ju =  cellfun(@(Jz,mub,mul)(Jz+(mub-mul)'),objPartials.Ju,muH.ub.u',muH.lb.u','UniformOutput',false);
    if withAlgs
        gbar.Jv = cellfun(@(Jz,mub,mul)(Jz+(mub-mul)'),objPartials.Jv,muH.ub.v',muH.lb.v','UniformOutput',false);
    end
    

    
    
    
    if relax || k == 1

        if  k > opt.multiplierFree
            gbarLambda.Jx = gbar.Jx;
            gbarLambda.Ju = gbar.Ju;
            if withAlgs
                gbarLambda.Jv = gbar.Jv;
            end
            [~,~,~,~,lambdaX,lambdaV]= simulateSystemZ(u,x,v,ss,[],'simVars',simVars,'JacTar',gbarLambda,'withAlgs',withAlgs);

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
            
        else
            normInfLambda = -inf;

        end        
        
        
        if xi < 1
            % multiplier free approximations
            [gbarR,errorSum,crossProduct] = multiplierFreeApproxs(gbar,ax,av,xd,vd,w,du,xi,withAlgs);
            % calculate equality constraints penalty
            [rho,errorSumB,dualApproxB] = equalityConsPenalty(gbarR,errorSum,crossProduct,rho,rhoHat,errorSumB,dualApproxB,normInfLambda);
        else
            warning('xi == 1. The problem may be infeasible to solve');
        end
    end
    
    %% Merit function definition
    merit = @(f,dE,varargin) l1merit(f,dE,rho,varargin{:});
    % line function

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
        'xi',xi,...
        'withAlgs',withAlgs,...
        varargin{:});
   
    
    % do not perform a watch-dog step on the very first iteration! 
    if k<=1
        skipWatchDog = true;
    else
        skipWatchDog = false;
    end
    
    % Line-search 
    [l,~,~,~,xfd,vars,simVars,relax,returnVars,wentBack,debugInfo] = watchdogLineSearch(phi,relax,...
        'tau',opt.tauL,'eta',opt.eta,'kmax',opt.lkMax,'debugPlot',opt.debugLS,'debug',opt.debug,...
        'simVars',simVars,'curvLS',opt.curvLS,'returnVars',returnVars,'skipWatchDog',skipWatchDog,'maxStep',maxStep,'k',k);
    
    
    if relax == false && (debugInfo{2}.eqNorm1 > debugInfo{1}.eqNorm1)  %% Watchdog step activated, should we perform SOC?
        
        
        % build the new problem!
        xsSOC =  bringVariables(vars.xs,jobSchedule);
        xSOC =  bringVariables(vars.x,jobSchedule);
        if withAlgs
            vsSOC = bringVariables(vars.vs,jobSchedule);
            vSOC = bringVariables(vars.v,jobSchedule);
        end
        
        xdSoc = cellfun(@(vxsi,vxi,xdi)vxsi-vxi+(1-xi)*xdi,xsSOC,xSOC,xd,'UniformOutput',false);
        if withAlgs
            vdSoc = cellfun(@(vvsi,vvi,vdi)vvsi-vvi+(1-xi)*vdi,vsSOC,vSOC,vd,'UniformOutput',false);
        else
            vdSoc = [];
        end
                
        [~,~,~,~,axSOC,~,avSOC,~] = condensing(x,u,v,ss,'simVars',simVars,'computeCorrection',true,'computeNullSpace',false,'xd',xdSoc,'vd',vdSoc,'withAlgs',withAlgs);

        if opt.computeCrossTerm
            [wSOC,stepYSOC] = computeCrossTerm(x,u,v,axSOC,avSOC,gbarZ,ss,obj,mudx,mudu,mudv,opt.lbxH,opt.lbvH,opt.ubxH,opt.ubvH,withAlgs,'xs',xs.client,'vs',vs.client);
        else
            wSOC = cellfun(@(xx)zeros([size(xx,2),size(xx,1)]),u','UniformOutput',false);
            stepYSOC = 0; 
        end
   
       
        if opt.condense
            upActiveSOC = opt.upActive;
            lowActiveSOC = opt.lowActive;
            Aact1 = [];
        end
        
        % Solve the QP to obtain the step on the nullspace.
        [ duSOC,dxSOC,dvSOC,xiSOC,lowActiveSOC,upActiveSOC,muHSOC,violationSOC,qpVAlSOC,dxNSOC,dvNSOC,~,QPITSOC] = qpStep(M,gZ,wSOC,...
            ldu,udu,...
            Aact1,predictor,constraintBuilder,...                      % TODO: reuse Aact from the convergence of the first QP, we just need to interface it!
            axSOC,Ax,ldx,udx,...
            avSOC,Av,ldv,udv,...
            'lowActive',lowActiveSOC,'upActive',upActiveSOC,...
            'ci',ss.ci,...
            'qpDebug',opt.qpDebug,'it',k,'withAlgs',withAlgs,...
            'feasTol',opt.qpFeasTol,'condense',opt.condense);

        QPIT = QPIT + QPITSOC;
        
        % debug check-point, check if the file is present
        if opt.debug
            fid = fopen('deleteMe2Break.txt','r');
            if fid == -1
                fid = fopen('deleteMe2Break.txt','w');fclose(fid);
                keyboard;
            else
                fclose(fid);
            end
        end
        
        if violationSOC.x > masterTol || (withAlgs && (violationSOC.v > masterTol))
            warning('QP solver too inaccurate, check the scaling and tolerance settings');
        end
        
        % Honor hard bounds in every step. Cut step if necessary
        [maxStep,duSOC] = maximumStepLength(u,duSOC,opt.lbu,opt.ubu,'tol',opt.qpFeasTol);
        
        [maxStepx,dxSOC] = maximumStepLength(x,dxSOC,opt.lbx,opt.ubx,'tol',violationSOC.x);
        maxStep = min(maxStep,maxStepx);
        if withAlgs
            [maxStepv,dvSOC] =maximumStepLength(v,dvSOC,opt.lbv,opt.ubv,'tol',violationSOC.v);
            maxStep = min(maxStep,maxStepv);
        end
        

        dxW = distributeVariables(dxSOC,jobSchedule);
        dvW = distributeVariables(dvSOC,jobSchedule);
        
        %trystep
        % TODO: implement function without calculating gradients!
        [ fSOC,dfSOC,varsSOC,simVarsSOC,debugInfoSOC ] = lineFunctionWrapper(maxStep,...
            xW,...
            vW,...
            u,...
            dxW,...
            dvW,...
            duSOC,...
            simFunc,obj,merit,jobSchedule,'gradients',true,'plotFunc',opt.plotFunc,'plot',opt.plot,...
            'withAlgs',withAlgs,...
            'xi',xi,...
            'debug',opt.debug);
        
        
        
        
        %Try full Step
        
        xfd = [xfd;maxStep fSOC dfSOC];
        
        armijoF = @(lT,fT)  (fT - (xfd(1,2) + opt.eta*xfd(1,3)*lT));
        armijoOk = @(lT,fT) (armijoF(lT,fT) <= 0);
        
        
        debugInfoSOC.armijoVal = armijoF(maxStep,fSOC);
        debugInfo = [debugInfo;debugInfoSOC];
              
        
        if armijoOk(maxStep,fSOC)  %% accept this step!
            ax = axSOC;
            av = avSOC;
            du = duSOC;
            dx = dxSOC;
            dv = dvSOC;
            xi = xiSOC;
            opt.lowActive = lowActiveSOC;
            opt.upActive = upActiveSOC;
            muH = muHSOC;
            violation = violationSOC;
            qpVAl = qpVAlSOC;
            dxN = dxNSOC;
            dvN = dvNSOC;
            w = wSOC;
            l=maxStep;
            vars = varsSOC;
            simVars = simVarsSOC;
            relax = true;
            returnVars = [];
            wentBack = false;
                        
            gbar.Jx =  cellfun(@(Jz,mub,mul)(Jz+(mub-mul)'),objPartials.Jx,muH.ub.x',muH.lb.x','UniformOutput',false);
            gbar.Ju =  cellfun(@(Jz,mub,mul)(Jz+(mub-mul)'),objPartials.Ju,muH.ub.u',muH.lb.u','UniformOutput',false);
            if withAlgs
                gbar.Jv = cellfun(@(Jz,mub,mul)(Jz+(mub-mul)'),objPartials.Jv,muH.ub.v',muH.lb.v','UniformOutput',false);
            end
                       
            
            debugWatchdog( k,'C',xfd(end,1),xfd(end,2),xfd(end,3),debugInfo(end));
        else
            debugWatchdog( k,'X',xfd(end,1),xfd(end,2),xfd(end,3),debugInfo(end));
        end
        
        
    
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
    
    % Restore previous lagrange multiplier estimate, if the Watch-dog
    % returned from a previous estimate
    if wentBack 
        mudx = muReturnX;
        if withAlgs
            mudv = muReturnV;
        end
        mudu = muReturnU;
        muH = muHReturn;
        gbarZm = gbarZmReturn;
        gbarZ = gbarZReturn;
   
    else
      
        % calculate the lagrangian with the updated values of mu, this will
        % help to perform the BFGS update
        if opt.condense
            if withAlgs
                gbarZm = vectorTimesZ(gbar.Jx,gbar.Ju,gbar.Jv,Ax,Av,ss.ci );
            else
                gbarZm = vectorTimesZ(gbar.Jx,gbar.Ju,[],Ax,[],ss.ci );
            end
        else
            [~,gbarZm,~,~,~,~] = simulateSystemZ(u,x,v,ss,[],'simVars',simVars,'JacTar',gbar,'withAlgs',withAlgs);
        end 
    end
    % Save Lagrange multipliers to restore if necessary
    if ~isempty(returnVars) && ~relax
        muReturnX = mudx;
        if withAlgs
            muReturnV = mudv;
        end
        muReturnU = mudu;
        muHReturn = muH;
        gbarZmReturn = gbarZm;
        gbarZReturn = gbarZ;
    else
        muReturnX = [];
        muReturnU = [];
        muReturnV = [];
        muHReturn = [];
        gbarZmReturn = [];
        gbarZReturn = [];
    end

    mudx = cellfun(@(x1,x2,x3)(1-l)*x1+l*(x2-x3),mudx,muH.ub.x,muH.lb.x,'UniformOutput',false);
    mudu = cellfun(@(x1,x2,x3)(1-l)*x1+l*(x2-x3),mudu,muH.ub.u,muH.lb.u,'UniformOutput',false);
    if withAlgs
        mudv = cellfun(@(x1,x2,x3)(1-l)*x1+l*(x2-x3),mudv,muH.ub.v,muH.lb.v,'UniformOutput',false);
    end
	gbarZm = cellfun(@(Z,Zm)(1-l)*Z+l*Zm,gbarZ,gbarZm,'UniformOutput',false);

    
    if l == 1
        wbar = w;
    else
        wbar = cellfun(@(wi)l*wi,w,'UniformOutput',false);
    end
    
    
    if opt.debug
        printLogLine(k,...
            {'|(g+nu)Z|','|c|','|Ypy|','|Zpz|','xi','|gZ|','|w|','pz''w','stepY','l','|s|','|y|','s''y','cond(B)'},...
            {...
            sqrt(sum(cellfun(@(x)dot(x,x),gbarZm))),...
            sqrt(sum(cellfun(@(x)dot(x,x),[xd;vd]))),...
            sqrt(sum(cellfun(@(x)dot(x,x),[ax;av]))),...
            sqrt(sum(cellfun(@(x)dot(x,x),[dxN;dvN;du]))),...
            xi,...
            sqrt(sum(cellfun(@(x)dot(x,x),gZ))),...
            sqrt(sum(cellfun(@(x)dot(x,x),w))),...
            sum(cellfun(@mtimes,w,du')),...
            stepY,...
            l,...
            norm(s),...
            norm(y),...
            sTy,...
            cond(M)...
            }...
            );
    end
    
    % save last value of u for the BFGS update
    uBV = cell2mat(u);
    
    % return the new iterate returned after line-search.
    x = bringVariables(vars.x,jobSchedule);  
    xs.worker = vars.xs;
    if isfield(xs,'client')
        xs = rmfield(xs,'client');
    end
    if withAlgs
        v =  bringVariables(vars.v,jobSchedule);
        vs.worker = vars.vs;
        if isfield(vs,'client')
            vs = rmfield(vs,'client');
        end
    end
    u = vars.u;

	[~,u]  = checkBounds( opt.lbu,u,opt.ubu,'chopp',true,'verbose',opt.debug);
    [~,x]  = checkBounds( opt.lbx,x,opt.ubx,'chopp',true,'verbose',opt.debug);
    if withAlgs
        [~,v]  = checkBounds( opt.lbv,v,opt.ubv,'chopp',true,'verbose',opt.debug);
    end    
    uV = cell2mat(u);

    usliced = [];
    
    % Save the current iteration to a file, for debug purposes.
    if opt.saveIt
        save itVars x u v xd vd rho M;
    end
    if ~isempty(opt.controlWriter)
        opt.controlWriter(u,k);
    end
    
    
    % print main debug
    if  opt.debug
        if mod(k,10) == 1
            header = true;
        else
            header = false;
        end
        tMax = 0;
        
        violationH = violation.x;
		if withAlgs
			violationH = max(violationH,violation.v);
		end
        dispFunc(k,norm(cell2mat(gbarZm)),violationH,normdu,rho,tMax,xfd,cond(M),relax,debugInfo,header,QPIT );
    end
    
    if l == 0  %line search couldn't make a sufficient decrease
        warning('lineSearch determined 0 step length');
        break;
    end
    
    
end


% recover previous variables if you performed a watch-dog step in the last iteration
if ~converged &&  ~relax
        x = bringVariables(returnVars.vars0.x,jobSchedule);
    u = returnVars.vars0.u;
    if withAlgs
            v = bringVariables(returnVars.vars0.v,jobSchedule);
    end
    simVars = returnVars.simVars0;
    [xs,vs,~,~,simVars] = simulateSystem(x,u,ss,'guessV',v,'simVars',simVars,'withAlgs',withAlgs);
    f = obj(xs,u,v,'gradients',false);
        xsF = bringVariables(xs,jobSchedule);
        xd = cellfun(@(x1,x2)x1-x2,xsF,x,'UniformOutput',false);
end
    
end


