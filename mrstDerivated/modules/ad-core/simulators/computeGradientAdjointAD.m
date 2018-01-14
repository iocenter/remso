function varargout = computeGradientAdjointAD(state0, states, model, schedule, getObjective, varargin)
%Compute gradients using an adjoint/backward simulation that is linear in each step
%
% SYNOPSIS:
%   grad = computeGradientAdjointAD(state0, states, model, schedule, getObjective
%
% DESCRIPTION:
%   For a given schedule, compute gradients with regards to well controls
%   by perturbing all controls ever so slightly and re-running the
%   simulation.
%  
%   As the cost of this routine grows is approximately
%
%        (# wells)x(# ctrl step) x cost of schedule
%
%   it can be potentially extremely expensive. It is better to use the
%   'computeGradientAdjointAD' routine for most practical purposes. This
%   routine is primarily designed for validation of said routine.
%
% REQUIRED PARAMETERS:
%
%   state0       - Physical model state at t = 0
%
%   states       - All previous states. Must support the syntax 
%                  state = states{i}. If the problem is too large to fit in
%                  memory, it can use ResultHandler class to retrieve files
%                  from the disk.
%                  
%   model        - Subclass of PhysicalModel class such as
%                 ThreePhaseBlackOilModel that models the physical effects
%                 we want to study.
%
%   schedule     - Schedule suitable for simulateScheduleAD.
%
%   getObjective - Function handle for getting objective function value 
%                  for a given timestep with derivatives. Format: @(tstep)
%
% OPTIONAL PARAMETERS (supplied in 'key'/value pairs ('pn'/pv ...)):
%   
% 
%
% RETURNS:
%
%  'Verbose'        - Indicate if extra output is to be printed such as
%                     detailed convergence reports and so on. 
% 
%  'scaling'        - Struct with fields 'rate' and 'pressure' used to
%                     scale the relevant control equations, if the model 
%                     supports it.
%
%  'LinearSolver'   - Subclass of 'LinearSolverAD' suitable for solving the
%                     adjoint systems.
%
% SEE ALSO:
%   computeGradientPerturbationAD, simulateScheduleAD

%{
Copyright 2009-2016 SINTEF ICT, Applied Mathematics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}

%{
Changes from https://github.com/iocenter/remso/ 4f8fa54a92c14b117445ad54e9e5dd3a0e47a7f5

Handle option of gradient*matrix 
Sensitivity with respect to initial conditions

%}
    opt = struct('ControlVariables', {{'well'}}, ...
                 'LinearSolver',     [],...
                 'xRightSeeds',[],...
                 'uRightSeeds',[]);
    opt = merge_options(opt, varargin{:});
    
    getState = @(i) getStateFromInput(schedule, states, state0, i);
    
    if isempty(opt.LinearSolver)
        linsolve = BackslashSolverAD();
    else
        linsolve = opt.LinearSolver;
    end
    
    if iscell(opt.ControlVariables) && strcmp(opt.ControlVariables{1},'well')
        ncv = 1;
    else
        error('Not supported, revise the code!')
    end
    
    nstep = numel(schedule.step.val);
    grad = [];
    gradstep = cell(nstep, ncv);
    nt = nstep;
    for step = nt:-1:1
        [dg, grad, report,gradstep{step}] = model.solveAdjoint(linsolve, getState, ...
                                         getObjective, schedule, grad, step);
    

    end
    
    
    % Sum up to the control steps
    nc = numel(schedule.control);
    gradients = cell(ncv, nc);
    for k = 1:nc
        ck = schedule.step.control == k;
        for j = 1:ncv
            tmp = vertcat(gradstep{ck, j});
            gradients{j, k} = full(sum(cat(3,tmp{:}), 3));
        end
    end
    if (size(opt.uRightSeeds,1)) ~= 0 || (size(opt.xRightSeeds,1) ~= 0)
        % conventions conventions...  should jac have variables on the columns or lines?
        gradients = cellfun(@(xit)xit',gradients,'UniformOutput',false);
    end
    
    if size(opt.uRightSeeds,1) ~= 0
        gradients = cell2mat(gradients)*cell2mat(opt.uRightSeeds);
    end
    if size(opt.xRightSeeds,1) ~= 0
        % observe that we could have given xRightSeeds to initialize the
        % function!
        [result0] = solveAdjointStete0Sens(model, getState,schedule, grad);
        if iscell(gradients)
            gradients = result0;
        else
            gradients = gradients+result0(:,1:size(opt.xRightSeeds,1))*opt.xRightSeeds;
        end
    end
    varargout{1} = gradients;
    if nargout >=2
        result0 = solveAdjointStete0Sens(model, getState,schedule, grad);
        nx = model.G.cells.num*sum(model.getActivePhases());
        varargout{2} = full(result0(:,1:nx)');
    end
end

function state = getStateFromInput(schedule, states, state0, i)
    if i == 0
        state = state0;
    elseif i > numel(schedule.step.val)
        state = [];
    else
        state = states{i};
    end
end
