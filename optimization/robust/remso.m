function [u,x,v,s,f,M,simVars] = remso(u,ss,obj,varargin)
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
    'qpDebug',true,...
    'lowActive',[],'upActive',[],...
    'simVars',[],'debug',true,'plot',false,'saveIt',false,...
    'controlWriter',[],...
    'multiplierFree',inf,...
    'allowDamp',true,...
    'qpFeasTol',1e-6,...
    'etaRisk',0.9,...
    'computeCrossTerm',true,...
    'condense',false,'testQP',false);

opt = merge_options(opt, varargin{:});


masterTol = min([opt.tol,opt.tolU,opt.tolX,opt.tolV]);

%The qpFeasTol must be tighter than tol, tolX, tolV, and tolU'
if opt.qpFeasTol > masterTol
    opt.qpFeasTol = masterTol;
end


% number of variables
nR = numel(ss);
xDims = cellfun(@(ssi)numel(ssi.state),ss);  %% it is assumed that the models can have different number of states
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

% contrasting the other implementations, use a single lbx for all the
% prediction horizon.
if isempty(opt.lbx)
    opt.lbx = arrayfun(@(nxi)-inf(nxi,1),xDims,'UniformOutput',false);
end
if isempty(opt.ubx)
    opt.ubx = arrayfun(@(nxi)-inf(nxi,1),uDims,'UniformOutput',false);
end
%% initial simulation profile
if isempty(opt.simVars)
    simVars = cell(nR,1);  % each cell inside most
else
    simVars = opt.simVars;
end

%% Process initial MS simulation guess, if not given, get it by forward simulation
simulateSS = false;
if ~isempty(opt.x)
    %  Initial guess for prediction given by the user
    x = opt.x;
    xs = opt.x;
else
    % Initial guess not provided, take from a simulation in the gradient
    % routine
    simulateSS = true;
    x = [];
    xs = [];
end

if isempty(opt.v)
    vs = [];
else
    vs = opt.v;
end

if simulateSS
    [~,~,~,simVars,xs,vs,s2,usliced] = simulateSystemSS_R(u,ss,[],'guessX',xs,'guessV',vs,'simVars',simVars);
    x = xs;
    v = vs;
    s = s2;
else
    [xs,vs,s2,~,~,simVars,usliced] = simulateSystem_R(x,u,v,ss,'gradients',false,'guessX',xs,'guessV',vs,'simVars',simVars);
    v = vs;
    s = s2;
end

xDims = cellfun(@(z)cellfun(@numel,z),x,'UniformOutput',false);
vDims = cellfun(@(z)cellfun(@numel,z),v,'UniformOutput',false);
assert(sum([0;cell2mat(vDims)])>0,'The robust optimization formulation must contain algebraic variables to be well defined');
withAlgs = true;

if isempty(opt.lbv)
    opt.lbv = cellfun(@(nz)arrayfun(@(nzk)-inf(nzk,1),nz,'UniformOutput',false),vDims,'UniformOutput',false);
end
if isempty(opt.ubv)
    opt.ubv = cellfun(@(nz)arrayfun(@(nzk) inf(nzk,1),nz,'UniformOutput',false),vDims,'UniformOutput',false);
end
if isempty(opt.ubs)
    opt.ubs = inf(size(s));
end
if isempty(opt.lbs)
    opt.lbs = -inf(size(s));
end

[~,x]  = cellfun(@(xr)    checkBounds( opt.lbx,xr,opt.ubx,'chopp',true,'verbose',opt.debug),        x,        'UniformOutput',false);
[~,v]  = cellfun(@(l,vr,u)checkBounds( l      ,vr,u,      'chopp',true,'verbose',opt.debug),opt.lbv,v,opt.ubv,'UniformOutput',false);
[~,s]  =                  checkBounds(opt.lbs ,s, opt.ubs,'chopp',true,'verbose',opt.debug);


Aact1= [];
predictor = [];
constraintBuilder = [];

% Multiple shooting simulation function
simFunc = @(xk,uk,vk,varargin) simulateSystem_R(xk,uk,vk,ss,'eta',opt.etaRisk,varargin{:});


%% Define empty active sets if they are not given
if isempty(opt.lowActive)
    opt.lowActive.x = cellfun(@(xr)cellfun(@(xrk)false(size(xrk)),xr,'UniformOutput',false),x,'UniformOutput',false);
    opt.lowActive.v = cellfun(@(vr)cellfun(@(vrk)false(size(vrk)),vr,'UniformOutput',false),v,'UniformOutput',false);
    opt.lowActive.s = {false(size(s))};
end
if isempty(opt.upActive)
    opt.upActive.x = cellfun(@(xr)cellfun(@(xrk)false(size(xrk)),xr,'UniformOutput',false),x,'UniformOutput',false);
    opt.upActive.v = cellfun(@(vr)cellfun(@(vrk)false(size(vrk)),vr,'UniformOutput',false),v,'UniformOutput',false);
    opt.upActive.s = {false(size(s))};
end



%% lagrange multipliers estimate initilization
mudx=  cellfun(@(xr)cellfun(@(xrk)sparse(1,numel(xrk)),xr','UniformOutput',false),x','UniformOutput',false);
mudv = cellfun(@(vr)cellfun(@(vrk)sparse(1,numel(vrk)),vr','UniformOutput',false),v','UniformOutput',false);
mudu = mat2cell(sparse(1,sum(uDims)),1,uDims);
muds = mat2cell(sparse(1,numel(s)),1,numel(s));

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
su = zeros(sum(uDims),1);
sTy = 0;

S = [];
Y = [];



%% Line-search parameters
rho = 1e-5;
rhoHat = 1e-5;
returnVars = [];
relax = false;   % to avoid the hessian update and perform a fine line-search
errorSumB = [];
dualApproxB = [];


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
    [xs,vs,s2,xd,vd,sd,ax,Ax,av,Av,as,As]  = condensing_R(x,u,v,s,ss,'simVars',simVars,'computeCorrection',true,'eta',opt.etaRisk,'computeNullSpace',opt.condense);
    
    [f,objPartials] = obj(s,u,'gradients',true);
    
    gbar.Ju = cellfun(@plus,objPartials.Ju,mudu,'UniformOutput',false);
    gbar.Js = objPartials.Js+cell2mat(muds);
    gbar.Jx = mudx;
    gbar.Jv = mudv;
    
    
    if opt.condense
        gZ = mat2cell(objPartials.Js*cell2mat(As) + cell2mat(objPartials.Ju),1,uDims);
        
        
        gbarZ = cellfun(@(gx,gv,gz,gu)gx+gv+gz+gu,...
            calcgbarZ(gbar.Jx,Ax,ss),...
            calcgbarZ(gbar.Jv,Av,ss),...
            mat2cell(gbar.Js*cell2mat(As),1,uDims),...
            gbar.Ju,...
            'UniformOutput',false);
        
    else
        objPartials.Jx = [];
        objPartials.Jv = [];
        
        
        [ sensitivities ] = generateSimulationSentivity(u,x,v,ss,simVars,[objPartials;gbar],xDims,vDims,uDims,opt.lowActive,opt.upActive );
        
        
        gZ = sensitivities{1};
        gbarZ = sensitivities{2};
        Aact1 = sensitivities{3};
        
        
        lowActiveSOC = opt.lowActive;
        upActiveSOC = opt.upActive;
        % if SOC is executed, start it with Aact1. TODO: use the jacobians
        % from the last QP, the later provides more information
        
        predictor = @(du) linearPredictor(du,x,u,v,s,ss,simVars);
        constraintBuilder = @(activeSet) generateSimulationSentivity(u,x,v,ss,simVars,[],xDims,vDims,uDims,activeSet);
        
        lagFunc = @(J)simulateSystemZ_R(u,x,v,ss,J,'simVars',simVars,'eta',opt.etaRisk);
        
        
    end
    
    
    
    
    
    
    
    if opt.computeCrossTerm
        
        % TODO: after finished check all input and outputs, in particular
        % uSliced!
        
        
        % Honor hard bounds in every step. Cut step if necessary
        [w,stepY] = computeCrossTerm(x,u,v,s,ax,av,as,gbarZ,ss,obj,mudx,mudu,mudv,muds,opt.lbx,opt.lbv,opt.lbs,opt.ubx,opt.ubv,opt.ubs,'xs',xs,'vs',vs,'s2',s2);
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
    
    
    %% Compute search direction  && lagrange multipliers
    
    % Compute bounds for the linearized problem
    udu =  cellfun(@(w,e)(w-e),opt.ubu,u,'UniformOutput',false);
    ldu =  cellfun(@(w,e)(w-e),opt.lbu,u,'UniformOutput',false);
    
    udx =  cellfun(@(er)cellfun(@(w,e)(w-e),opt.ubx,er,'UniformOutput',false),x,'UniformOutput',false);
    ldx =  cellfun(@(er)cellfun(@(w,e)(w-e),opt.lbx,er,'UniformOutput',false),x,'UniformOutput',false);
    udv =  cellfun(@(wr,er)cellfun(@(w,e)(w-e),wr,er,'UniformOutput',false),opt.ubv,v,'UniformOutput',false);
    ldv =  cellfun(@(wr,er)cellfun(@(w,e)(w-e),wr,er,'UniformOutput',false),opt.lbv,v,'UniformOutput',false);
    uds =  opt.ubs-s;
    lds =  opt.lbs-s;
    
    % Solve the QP to obtain the step on the nullspace.
    [ du,dx,dv,ds,xi,opt.lowActive,opt.upActive,muH,violation,qpVAl,dxN,dvN,dsN,slack,QPIT] = qpStep_R(M,gZ,w,...
        ldu,udu,...
        Aact1,predictor,constraintBuilder,...
        ax,Ax,ldx,udx,...
        av,Av,ldv,udv,...
        as,As,lds,uds,...
        'lowActive',opt.lowActive,'upActive',opt.upActive,...
        'ss',ss,...
        'qpDebug',opt.qpDebug,'it',k,...
        'feasTol',opt.qpFeasTol,'condense',opt.condense,'lagFunc',lagFunc,'testQP',opt.testQP);
    
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
    
    if (violation.x > masterTol) || (violation.v > masterTol) || ((violation.s > masterTol))
        warning('QP solver too inaccurate, check the scaling and tolerance settings');
    end
    
    
    % Honor hard bounds in every step. Cut step if necessary, use the QP
    % tolerance setting to do so
    [maxStep,du] = maximumStepLength(u,du,opt.lbu,opt.ubu,'tol',opt.qpFeasTol);
    
    [maxStepx,dx] = cellfun(@(zi,dz)maximumStepLength(zi,dz,opt.lbx,opt.ubx,'tol',violation.x),x,dx,'UniformOutput',false);
    maxStep = min(maxStep,max(cell2mat(maxStepx)));
    [maxStepv,dv] = cellfun(@(zi,dz,lb,ub)maximumStepLength(zi,dz,lb,ub,'tol',violation.v),v,dv,opt.lbv,opt.ubv,'UniformOutput',false);
    maxStep = min(maxStep,max(cell2mat(maxStepv)));
    [maxSteps,ds] = maximumStepLength({s},{ds},{opt.lbs},{opt.ubs},'tol',violation.s);
    ds = cell2mat(ds);
    maxStep = min(maxStep,max(maxSteps));
    
    
    
    
    %% Convergence test
    % I choose the infinity norm, because this is easier to relate to the
    % physical variables
    normdu = norm(cellfun(@(z)norm(z,'inf'),du),'inf');
    normax = norm(cellfun(@(z)norm(cell2mat(z),'inf'),ax),'inf');
    normav = norm(cellfun(@(z)norm(cell2mat(z),'inf'),av),'inf');
    normas = norm(as,'inf');
    
    if normdu < opt.tolU && normax < opt.tolX && normav < opt.tolV && normas < opt.tolS  && normdu < opt.tol && normax < opt.tol && normav < opt.tol &&  normas < opt.tol && relax
        converged = true;
        break;
    end
    
    %% Preparing for line-search
    
    
    
    % gbar = g+nu
    gbar.Ju = cellfun(@plus,objPartials.Ju,minusS(muH.ub.u,muH.lb.u),'UniformOutput',false);
    gbar.Js = objPartials.Js+cell2mat(minusS(muH.ub.s,muH.lb.s));
    gbar.Jx = minusC(muH.ub.x,muH.lb.x);
    gbar.Jv = minusC(muH.ub.v,muH.lb.v);
    
    if relax || k == 1
        
        if  k > opt.multiplierFree
            gbarLambda.Jx = gbar.Jx;
            gbarLambda.Ju = gbar.Ju;
            gbarLambda.Jv = gbar.Jv;
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
        
        
        if xi ~= 1
            % multiplier free approximations
            [gbarR,errorSum,crossProduct] = multiplierFreeApproxs(gbar,ax,av,as,xd,vd,sd,w,du,xi);
            % calculate equality constraints penalty
            [rho,errorSumB,dualApproxB] = equalityConsPenalty(gbarR,errorSum,crossProduct,rho,rhoHat,errorSumB,dualApproxB,normInfLambda);
        else
            warning('xi == 1. The problem may be infeasible to solve');
        end
    end
    
    %% Merit function definition
    merit = @(f,dE,dS,varargin) l1merit(f,dE,dS,rho,varargin{:});
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
        simFunc,obj,merit,'gradients',true,'plotFunc',opt.plotFunc,'plot',opt.plot,...
        'debug',opt.debug,...
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
        'tau',opt.tauL,'eta',opt.eta,'kmax',opt.lkMax,'debugPlot',opt.debugLS,'debug',opt.debug,...
        'simVars',simVars,'curvLS',opt.curvLS,'returnVars',returnVars,'skipWatchDog',skipWatchDog,'maxStep',maxStep,'k',k);
    
    
    if relax == false && (debugInfo{2}.eqNorm1 > debugInfo{1}.eqNorm1)  %% Watchdog step activated, should we perform SOC?
        
        
        % build the new problem!
        xdSOC =  zdSOC(vars.xs,vars.x,xd,xi);
        vdSOC =  zdSOC(vars.vs,vars.v,vd,xi);
        sdSOC = zdSOCS({vars.s2},{vars.s},{sd},xi);
        sdSOC = cell2mat(sdSOC);
        
        
        [~,~,~,~,~,~,axSOC,~,avSOC,~,asSOC,~] = condensing_R(x,u,v,s,ss,...
            'simVars',simVars,...
            'computeCorrection',true,...
            'computeNullSpace',false,...
            'xd',xdSOC,'vd',vdSOC,'sd',sdSOC,...
            'eta',opt.etaRisk);
        
        
        if opt.computeCrossTerm
            [wSOC,stepYSOC] = computeCrossTerm(x,u,v,s,...
                axSOC,avSOC,asSOC,...
                gbarZ,ss,obj,...
                mudx,mudu,mudv,muds,...
                opt.lbx,opt.lbv,opt.lbs,opt.ubx,opt.ubv,opt.ubs,...
                'xs',xs,'vs',vs,'s2',s2);
        else
            stepYSOC = 0;
            wSOC = mat2cell(zeros(1,sum(uDims)),1,uDims);
        end
        
        
        if opt.condense
            upActiveSOC = opt.upActive;
            lowActiveSOC = opt.lowActive;
            Aact1 = [];
        end
        
        % Solve the QP to obtain the step on the nullspace.
        [ duSOC,dxSOC,dvSOC,dsSOC,xiSOC,lowActiveSOC,upActiveSOC,muHSOC,violationSOC,qpVAlSOC,dxNSOC,dvNSOC,dsNSOC,slack,QPITSOC] = qpStep_R(M,gZ,wSOC,...
            ldu,udu,...
            Aact1,predictor,constraintBuilder,...
            axSOC,Ax,ldx,udx,...
            avSOC,Av,ldv,udv,...
            asSOC,As,lds,uds,...
            'lowActive',lowActiveSOC,'upActive',upActiveSOC,...
            'ss',ss,...
            'qpDebug',opt.qpDebug,'it',k,...
            'feasTol',opt.qpFeasTol,'condense',opt.condense,'lagFunc',lagFunc,'testQP',opt.testQP);
        
        QPIT = QPIT+ QPITSOC;
        
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
        
        
        if (violationSOC.x > masterTol) || (violationSOC.v > masterTol) || ((violationSOC.s > masterTol))
            warning('QP solver too inaccurate, check the scaling and tolerance settings');
        end
        
        
        % Honor hard bounds in every step. Cut step if necessary, use the QP
        % tolerance setting to do so
        [maxStep,duSOC] = maximumStepLength(u,duSOC,opt.lbu,opt.ubu,'tol',opt.qpFeasTol);
        
        [maxStepx,dxSOC] = cellfun(@(zi,dz)maximumStepLength(zi,dz,opt.lbx,opt.ubx,'tol',violationSOC.x),x,dxSOC,'UniformOutput',false);
        maxStep = min(maxStep,max(cell2mat(maxStepx)));
        [maxStepv,dvSOC] = cellfun(@(zi,dz,lb,ub)maximumStepLength(zi,dz,lb,ub,'tol',violationSOC.v),v,dvSOC,opt.lbv,opt.ubv,'UniformOutput',false);
        maxStep = min(maxStep,max(cell2mat(maxStepv)));
        [maxSteps,dsSOC] = maximumStepLength({s},{dsSOC},{opt.lbs},{opt.ubs},'tol',violationSOC.s);
        dsSOC = cell2mat(dsSOC);
        maxStep = min(maxStep,max(maxSteps));
        
        
        [ fSOC,dfSOC,varsSOC,simVarsSOC,debugInfoSOC ] = lineFunctionWrapper(maxStep,...
            x,...
            v,...
            u,...
            s,...
            dxSOC,...
            dvSOC,...
            duSOC,...
            dsSOC,...
            simFunc,obj,merit,'gradients',true,'plotFunc',opt.plotFunc,'plot',opt.plot,...
            'debug',opt.debug,...
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
            opt.lowActive = lowActiveSOC;
            opt.upActive = upActiveSOC;
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
            gbar.Ju = cellfun(@plus,objPartials.Ju,minusS(muH.ub.u,muH.lb.u),'UniformOutput',false);
            gbar.Js = objPartials.Js+cell2mat(minusS(muH.ub.s,muH.lb.s));
            gbar.Jx = minusC(muH.ub.x,muH.lb.x);
            gbar.Jv = minusC(muH.ub.v,muH.lb.v);
            
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
            
            gbarZm = cellfun(@(gx,gv,gz,gu)gx+gv+gz+gu,...
                calcgbarZ(gbar.Jx,Ax,ss),...
                calcgbarZ(gbar.Jv,Av,ss),...
                mat2cell(gbar.Js*cell2mat(As),1,uDims),...
                gbar.Ju,...
                'UniformOutput',false);
            
        else
            gbarZm = simulateSystemZ_R(u,x,v,ss,gbar,'simVars',simVars,'eta',opt.etaRisk);
            
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
    mudx = estimateDuals(mudx,muH.ub.x,muH.lb.x,l);
    mudv = estimateDuals(mudv,muH.ub.v,muH.lb.v,l);
    mudu = estimateDualsS(mudu,muH.ub.u,muH.lb.u,l);
    muds = estimateDualsS(muds,muH.ub.s,muH.lb.s,l);
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
            sqrt(sum(dotSum(xd)+dotSum(vd)+sum(dot(sd,sd)))),...
            sqrt(sum(dotSum(ax)+dotSum(av)+sum(dot(as,as)))),...
            sqrt(sum(dotSum(dxN)+dotSum(dvN)+dotSumS({dsN})+dotSumS(du))),...
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
    
    % save last value of u for the BFGS update
    uBV = cell2mat(u);
    
    % return the new iterate returned after line-search.
    x = vars.x;
    xs = vars.xs;
    v =  vars.v;
    vs = vars.vs;
    s =  vars.s;
    s2 = vars.s2;
    
    u = vars.u;
    
    [~,x]  = cellfun(@(xr)    checkBounds( opt.lbx,xr,opt.ubx,'chopp',true,'verbose',opt.debug),        x,        'UniformOutput',false);
    [~,v]  = cellfun(@(l,vr,u)checkBounds( l      ,vr,u,      'chopp',true,'verbose',opt.debug),opt.lbv,v,opt.ubv,'UniformOutput',false);
    [~,s]  =                  checkBounds(opt.lbs ,s, opt.ubs,'chopp',true,'verbose',opt.debug);
    
    uV = cell2mat(u);
    
    usliced = [];
    
    % Save the current iteration to a file, for debug purposes.
    if opt.saveIt
        save itVars x u v xd vd vs rho M;
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
        violationH = max(violationH,violation.v);
        
        dispFunc(k,norm(cell2mat(gbarZm)),violationH,normdu,rho,tMax,xfd,cond(M),relax,debugInfo,header,QPIT );
    end
    
    if l == 0  %line search couldn't make a sufficient decrease
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

function me = estimateDuals(m,mU,mL,l)
f = @(mr,mUr,mLr)estimateDualsS(mr,mUr,mLr,l);
me = cellfun(f,m,mU,mL,'UniformOutput',false);
end

function me = estimateDualsS(m,mU,mL,l)
me = cellfun(@(x1,x2,x3)(1-l)*x1+l*(x2-x3),m,mU,mL,'UniformOutput',false);
end

function me = minusC(mU,mL)
f = @minusS;
me = cellfun(f,mU,mL,'UniformOutput',false);
end

function me = minusS(mU,mL)
me = cellfun(@(x1,x2)(x1-x2),mU,mL,'UniformOutput',false);
end

function me = calcgbarZ(J,A,ss)
f = @calcgbarZS;
vr = cellfun(f,J',A,ss,'UniformOutput',false);

me = catAndSum(vr);
end
function v = calcgbarZS(J,A,ss)
v = cellmtimesT( J,A,'lowerTriangular',true,'ci',ss.ci,'columnVector',false);
end

function out = catAndSum(M)

dims = cellfun(@(x)size(x,2),M{1});
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

out = mat2cell(out,size(out,1),dims);

end

function e = dotSum(z)
f = @dotSumS;
e = sum(cellfun(f,z));
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
