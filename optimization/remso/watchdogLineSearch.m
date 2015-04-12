function [l,fL,gL,neval,xfd,vars,simVars,nextRelax,returnVars,wentBack,debugInfo] = watchdogLineSearch(fun,relax,varargin)
% Watchdog line-search inspired in:
%
% [1] R. H. Byrd and J. Nocedal,
% “An analysis of reduced Hessian methods for constrained optimization,”
% Math. Program., vol. 49, no. 1–3, pp. 285–323, Nov. 1990.
%
%
% SYNOPSIS:
%
%  [l,fL,gL,neval,xfd,vars,simVars,nextRelax,returnVars,wentBack,debugInfo] = watchdogLineSearch(fun,relax)
%  [l,fL,gL,neval,xfd,vars,simVars,nextRelax,returnVars,wentBack,debugInfo] = watchdogLineSearch(fun,relax, 'pn', pv, ...)
%
% PARAMETERS:
%
%   fun - anonymous funtion fun: [0,1] -> R^{1}, the line function
%
%   relax - true if relaxed step is allowed.
%
%   'pn'/pv - List of 'key'/value pairs defining optional parameters. The
%             supported options are:
%
%   skipRelaxRatio - Parameter to skip a relaxed step
%
%   tau - parameter ralated to the gradient condition on the line-search
%
%   eta - parameter related to the armijo 'suficient-decrease' condition
%
%   kmax - maximum number of extra points to be evaluated
%
%   simVars - Simulation variables at the beguining of the line.
%
%   curvLS - if true then honor the curvature condition during line-search
%
%   returnVars - Variables to restore in case that the wathcdog fails
%
%   skipWatchDog - Perform normal line-search.
%
%   debug - true to collect debug information
%
%   maxStep - Maximum step length (default =1)
%
% RETURNS:
%
%   (l,fL,gL,) - Line value, function value and gradient obtained
%
%   neval - number of point evaluations during line-search
%
%   xfd - record of line value, function, and gradients during line-search
%
%  (vars,simVars) - inputs and outputs of fun at the line-search solution
%  that are returned to avoid recomputation
%
%   nextRelax - true if next step is relax, false otherwise.
%
%   returnVars - variables to restore after a failed watchdog step
%
%   wentBack - true if watchdog restored the previous step
%
%   debugInfo - debug information
%
%
% SEE ALSO:
%
%

opt = struct('skipRelaxRatio',5,'tau',0.1,'eta',1e-1,'kmax',4,'simVars',[],'curvLS',true,'returnVars',[],'skipWatchDog',false,'debugPlot',false,'debug',false,'maxStep',1,'k',0);
opt = merge_options(opt, varargin{:});

returnVars = [];
skipLineSearh = false;
wentBack = false;
nextRelax = true;
% Check descent condition


% evaluate the function at 0
[f0,g0,vars0,simVars0,debugInfoN]=  fun(0,'simVars',opt.simVars);

armijoF = @(lT,fT)  (fT - (f0 + opt.eta*g0*lT));
armijoOk = @(lT,fT) (armijoF(lT,fT) <= 0);



xfd = [0,f0,g0];
debugInfoN.armijoVal = armijoF(0,f0);
debugInfo = {debugInfoN};
neval = 1;

if  g0 >= 0

    
    l = 0;
    fL = f0;
    gL = g0;
    vars = vars0;
    simVars = simVars0;
    skipLineSearh = true;
    nextRelax = true;
   
    % TODO: What should be the behaviour if skipWatchDog == true during a
    % watchdog step? (i.e if there are variables to return to)
    if opt.skipWatchDog
		warning('check the penalty calculation procedure!');
        return;
    elseif relax == true  %%% there is no turning back-point!
        assert(isempty(opt.returnVars))
		warning('check the penalty calculation procedure!');
        return;
    else              %%% go back to the previous iteration
        fL = inf;
        gL = inf;
        wentBack = true;
    end
end



% if skipWatchDog, perform normal line-search and forget about everything
if opt.skipWatchDog
    %Try full Step
    l = opt.maxStep;        
    [fL ,gL,vars,simVars,debugInfoN] = fun(l);
    neval = neval +1;
    xfd = [xfd;l fL gL];
    
    
    
    debugInfoN.armijoVal = armijoF(l,fL);
    debugInfo = [debugInfo;debugInfoN];
    
    if ~armijoOk(l,fL)
        
        [l,fL,gL,vars,simVars,nevalN,xfdN,debugInfoN] = linesearch(fun,f0,g0,fL,gL,opt.eta,...
            'tau',opt.tau,...
            'kmax',opt.kmax,...
            'debug',opt.debugPlot,...
            'curvLS',opt.curvLS,...
            'maxStep',l);
        
        if l <= 0
            l = 0;
            fL = f0;
            gL = g0;
            vars = vars0;
            simVars= simVars0;
        end
        
        neval = neval + nevalN;
        xfd = [xfd;xfdN];
        debugInfo = [debugInfo;debugInfoN];
    end
    if opt.debug
        debugWatchdog( opt.k,'N',xfd(:,1),xfd(:,2),xfd(:,3),debugInfo);
    end
    return
end


if relax  % try to relax the line-search conditions
    
    %Try full Step
    l = opt.maxStep;        % try full
    [fL ,gL,vars,simVars,debugInfoN] = fun(l);
    neval = neval +1;
    xfd = [xfd;l fL gL];
    debugInfoN.armijoVal = armijoF(l,fL);
    debugInfo = [debugInfo;debugInfoN];
    
    if (fL <= f0 + opt.eta*g0*l);
        nextRelax = true;  %% relaxed step again!
        
        if opt.debug
            debugWatchdog( opt.k,'A',xfd(:,1),xfd(:,2),xfd(:,3),debugInfo);
        end
        
    elseif ((((f0-fL)/(g0*l) > opt.skipRelaxRatio))) || (opt.maxStep < 1)
        assert(~skipLineSearh,'Do not make a lineSerach with g0>=0');
        nextRelax = true;  %% too bad, skip watchdog step!
        
        [l,fL,gL,vars,simVars,nevalN,xfdN,debugInfoN] = linesearch(fun,f0,g0,fL,gL,opt.eta,...
            'tau',opt.tau,...
            'kmax',opt.kmax,...
            'debug',opt.debugPlot,...
            'curvLS',opt.curvLS,...
            'maxStep',l);
        if l <= 0
            l = 0;
            fL = f0;
            gL = g0;
            vars = vars0;
            simVars= simVars0;
        end
        
        neval = neval + nevalN;
        xfd = [xfd;xfdN];
        debugInfo = [debugInfo;debugInfoN];
        
        if opt.debug
            debugWatchdog( opt.k,'B',xfd(:,1),xfd(:,2),xfd(:,3),debugInfo);
        end
    else
        nextRelax = false;  %% Start watchdog - next step will be more careful :) !
        returnVars.fun = fun;
        returnVars.f0 = f0;
        returnVars.g0 = g0;
        returnVars.vars0 = vars0;
        returnVars.simVars0 = simVars0;
        returnVars.f1 = fL;
        returnVars.g1 = gL;
        returnVars.debugInfo = debugInfo;
        returnVars.xfd = xfd;
        
        if opt.debug
            debugWatchdog( opt.k,'W',xfd(:,1),xfd(:,2),xfd(:,3),debugInfo);
        end
    end
    
else % watchdog step - merit function decrease must be achieved!
    
    % now we need to recover the bad performance inthe previous step
    requiredDecrease =  opt.returnVars.f0 + opt.eta*opt.returnVars.g0;
    
    if ~skipLineSearh  %% positive search direction!, go back without trying line-search
        
        %Try full Step
        l = opt.maxStep;        
        [fL ,gL,vars,simVars,debugInfoN] = fun(l);
        neval = neval +1;
        xfd = [xfd;l fL gL];
        debugInfoN.armijoVal = armijoF(l,fL);
        debugInfo = [debugInfo;debugInfoN];
        
        % Sufficient decrease is required for this step and the step
        % before!
        if ~(armijoOk(l,fL) && (fL <= requiredDecrease))   
            
            % line-search requiring sufficient decrease for both steps
            [l,fL,gL,vars,simVars,nevalN,xfdN,debugInfoN] = linesearch(fun,f0,g0,fL,gL,opt.eta,...
                'tau',opt.tau,...
                'kmax',opt.kmax,...
                'debug',opt.debugPlot,...
                'curvLS',opt.curvLS,...
                'required',requiredDecrease,...
                'maxStep',l);
            
            if l <= 0
                l = 0;
                fL = f0;
                gL = g0;
                vars = vars0;
                simVars= simVars0;
            end
            
            neval = neval + nevalN;
            xfd = [xfd;xfdN];
            debugInfo = [debugInfo;debugInfoN];
            
        end
        
        
    end
    nextRelax = true;
    
    if (fL <= requiredDecrease) && ~skipLineSearh
       % GREAT SUCCESS!
        if opt.debug
            debugWatchdog( opt.k,'S',xfd(:,1),xfd(:,2),xfd(:,3),debugInfo);
        end
    else
        % Watchdog procedure failed, go back to the initial point
        % compute next step using normal line search from the first point
        % and finish the whatchdog
        
        [l,fL,gL,vars,simVars,nevalN,xfdN,debugInfoN] = linesearch(opt.returnVars.fun,...
            opt.returnVars.f0,...
            opt.returnVars.g0,...
            opt.returnVars.f1,...
            opt.returnVars.g1,...
            opt.eta,...
            'tau',opt.tau,...
            'kmax',opt.kmax,...
            'debug',opt.debugPlot,...
            'curvLS',opt.curvLS,...
            'maxStep',1);
        if l <= 0
            l = 0;
            fL = opt.returnVars.f0;
            gL = opt.returnVars.g0;
            vars = opt.returnVars.vars0;
            simVars= opt.returnVars.simVars0;
        end
        
        if opt.debug
            debugWatchdog( opt.k,'F',xfd(:,1),xfd(:,2),xfd(:,3),debugInfo);
            debugWatchdog( opt.k,'R',xfdN(:,1),xfdN(:,2),xfdN(:,3),debugInfoN,'header',false);
        end
        
        xfd = [xfd;opt.returnVars.xfd;xfdN];
        debugInfo = [debugInfo;opt.returnVars.debugInfo;debugInfoN];
        neval = neval+nevalN;
        wentBack = true;
        
    end
end



end

