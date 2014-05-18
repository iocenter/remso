function [state, meta,eqs] = stepOW(state0, state, meta, dt, G, W, system, fluid, varargin)
% Do a single step of a nonlinear solve for a Oil-Water system
% This function should in general not be called directly and is as such not
% documented with regards to input/output: See solvefiADI for an
% explanation of how the ad-fi solvers are implemented.

%{
#COPYRIGHT#
%}
%{
Modification by Codas:
* Inclusion of same changes suggested by Stein on the simulation solver
* include option 'iterative'
* linear solver for multiple rhs
* Output of the Jacobian
* Save information on convergence
* Take care of the 'iteration' option for compatibility
%}

% $Date: 2013-11-27 00:33:36 +0100 (on, 27 nov 2013) $
% $Revision: 12077 $

opt = struct('Verbose', mrstVerbose, 'temperature', false, 'minerals', false);
opt = merge_options(opt, varargin{:});
s = system.s;

% if ~isempty(system.podbasis)
%     solve = @(eqs) SolveEqsADIPOD(eqs, opt.podbasis);
% else

fprops = functions(system.getEquations);
if strcmp(fprops.function,'eqsfiOW') == 1
    [eqs, state] = system.getEquations(state0, state, dt, G, W, s, fluid, ...
        'temperature', opt.temperature, ...
        'minerals', opt.minerals, ...
        'iteration', meta.iteration);
    
else
    eqs = system.getEquations(state0, state, dt, G, W, s, fluid, ...
        'temperature', opt.temperature, ...
        'minerals', opt.minerals);
end

if system.nonlinear.cpr && isempty(system.podbasis)
    p  = mean(state0.pressure);
    bW = fluid.bW(p); bO = fluid.bO(p);
    sc = [1./bO, 1./bW];
    
    vargs = { 'ellipSolve', system.nonlinear.cprEllipticSolver, ...
        'cprType'   , system.nonlinear.cprType          , ...
        'relTol'    , system.nonlinear.cprRelTol        , ...
        'eqScale'   , sc                                , ...
        'iterative' , system.nonlinear.itLinearSolver};
    
    [dx, gmresits, solver_diverged] = cprGenericM(eqs, system, vargs{:});
    
else
    
    dx = SolveEqsADI(eqs, system.podbasis);
    solver_diverged = false;
    gmresits = [0 0];
end



% dx = solve(eqs);

% [state, nInc] = updateState(state, dx);


searchfail = true;
if system.nonlinear.linesearch
    getEqs = @(state) system.getEquations(state0, state, dt, G, W, s, fluid, 'resOnly', true,...
        'temperature', opt.temperature,...
        'minerals',opt.minerals);
    upState = @(dx) updateState(state, dx, opt);
    [state, dx, searchfail] = linesearchADI(state, dx, system, getEqs, upState, false);
end

% Update reservoir conditions once a delta has been found.
if searchfail
    [dx, meta] = stabilizeNewton(dx, meta, system);
    % If the line search failed, uncritically accept the first step and
    % pray the other measures (relaxation / dampening) handle the error.
    [state, nInc] = updateState(state, dx, opt);
end

[converged CNV MB] = getConvergence(state, eqs, fluid, system, dt);
[meta, residuals] = getResiduals(meta, eqs, system, solver_diverged);
meta.CNV = CNV;
meta.MB = MB;
if(opt.temperature || opt.minerals)
%{
    %residuals = cellfun(@(x) norm(x.val, 'inf'), eqs);
    %meta.converged = all(residuals < system.nonlinear.tol);
%}


    if opt.Verbose
        
        if meta.iteration == 1
            eqnnames = {'Oil', 'Water',  'qOs', 'qWs', 'pBHP'};
            if(opt.temperature)
                eqnnames{end+1} = 'T';
            end
            
            if(opt.minerals)
                eqnnames{end+1} = 'MI';
            end
            fprintf('%-9s', eqnnames{:})
            fprintf('\n');
        end
        fprintf('%8.2e ', residuals);
    end
else
        meta.converged = converged;
        meta.stopped = meta.iteration == system.nonlinear.maxIterations && ~converged;
        
    if opt.Verbose
%        residuals = cellfun(@(x) norm(x.val, 'inf'), eqs);
        eqnnames = {'Oil', 'Water',  'qOs', 'qWs', 'control'};
        printResidual(residuals, gmresits, eqnnames, meta.iteration, CNV, MB);
    end
end
end

%--------------------------------------------------------------------------

function [state, nInc] = updateState(state, dx, opt)
dsMax = .2;
dpMax = .3;
% if ~isempty(phi)
%     for i = 1:numel(dx)
%         dx{i} = phi.basis{i}*dx{i};
%     end
% end

dp = dx{1};
ds = dx{2};
nInc = max( norm(dp,'inf')/norm(state.pressure, 'inf'), ...
            norm(ds,'inf')/norm(state.s(:,1), 'inf') );
    
%maxch = norm(ds, 'inf');
%step = min(1, maxSatStep./maxch);

ds = sign(ds).*min(abs(ds), dsMax);
dp = sign(dp).*min(abs(dp), abs(dpMax.*state.pressure));

state.pressure = state.pressure + dp;
sw = state.s(:,1) + ds;
% Cap values
sw = min(sw, 1); sw = max(sw, 0);

state.s = [sw, 1-sw];

dqWs    = dx{3};
dqOs    = dx{4};
dpBHP   = dx{5};
var_num=5;
if(opt.temperature)
    var_num = var_num+1;
    state.T=state.T+step*dx{var_num};
    state.T = max(state.T,273);
    state.T = min(state.T,500);
end

if(opt.minerals)
    for i=1:size(state.I,2)
        var_num=var_num+1;
        state.I(:,i)=state.I(:,i)+dx{var_num};
        state.I(:,i)=max(state.I(:,i),0);
    end
    for i=1:size(state.M,2)
        var_num=var_num+1;
        state.M(:,i)=state.M(:,i)+dx{var_num};
        state.M(:,i)=max(state.M(:,i),0);
    end
end
if isfield(state.wellSol,'pressure')  %% Compatibility: Why did they change names!
    dpBHP = sign(dpBHP).*min(abs(dpBHP), abs(dpMax.*vertcat(state.wellSol.pressure)));
    for w = 1:numel(state.wellSol)
        state.wellSol(w).pressure      = state.wellSol(w).pressure + dpBHP(w);
        state.wellSol(w).qWs      = state.wellSol(w).qWs + dqWs(w);
        state.wellSol(w).qOs      = state.wellSol(w).qOs + dqOs(w);
        %    mw = min(1, max(0, state.wellSol(w).mixs(:,1) + dmixWs(w)));
        %    state.wellSol(w).mixs     = [mw, 1-mw];
    end
else
    dpBHP = sign(dpBHP).*min(abs(dpBHP), abs(dpMax.*vertcat(state.wellSol.bhp)));
    for w = 1:numel(state.wellSol)
        state.wellSol(w).bhp      = state.wellSol(w).bhp + dpBHP(w);
        state.wellSol(w).qWs      = state.wellSol(w).qWs + dqWs(w);
        state.wellSol(w).qOs      = state.wellSol(w).qOs + dqOs(w);
        %    mw = min(1, max(0, state.wellSol(w).mixs(:,1) + dmixWs(w)));
        %    state.wellSol(w).mixs     = [mw, 1-mw];
    end
end
end
