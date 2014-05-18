function [state, its, convergence,eqs] = ...
      solvefiADI(state0, dt, W, G, system, varargin)
%Solve single timestep of automatic differentiation system
%
% SYNOPSIS:
%   state = solvefiADI(state0, dt, W, G, system)
%   state = solvefiADI(state0, dt, W, G, system, ...)
%
% PARAMETERS:
%   state0 - Initial state as defined by initResSol
%
%   dt     - Size of timestep.
%
%   W      - Well configuration. Must have valid sign for all producers and
%            injectors.
%
%   G      - Grid. See grid_structure.
%
%   system - Valid system as defined by initADISystem. This function
%            defines all options for convergence, how to set up equations
%            and parameters for the non-linear solver.
%
%   ...    - Additional parameters passed on to 'system.stepFunction'.
%            OPTIONAL.
%
% RETURNS:
%   state  - Updated state after convergence or max iterations
%            being reached.
%
% SEE ALSO:
%   stepOW, stepBlackOil, stepOWPolymer, runScheduleADI
 
% COMMENTS:
%  
%   The fully implicit AD solvers are written to minimize duplication of
%   coding effort. Primarly, this consists of a splitting of general logic
%   for solving schedules and timesteps and the setup of equations and
%   their solution. In general, a system is created which defines the
%   phases present along with handles to the stepFunction for that specific
%   type of system.
%
%   When solving for a single timestep, the stepFunction is called until
%   the stepFunction reports convergence or stopping in the meta struct it
%   returns. 
%
%   The step functions can be completely arbitrary, but for the current
%   solvers the general structure are:
%
%    - solvefiADI calls stepFunction
%    - The stepFunction calls eqsfi<systemtype>, solves a single Newton
%      step to accuracies defined by the system input using some kind of
%      linear solver (Direct solver / some CPR preconditioner
%      configuration) and then updates the state.
%    - solvefiADI checks wells and returns the final state.
%
%   Depending on the simulation, solvefiADI may be exposed directly in a
%   for loop or hidden within a wrapper such as runScheduleADI.

%{
#COPYRIGHT#
%}
%{
Modification by Codas:
* Some changes by Stein
* initialGuess
* meta.iteration
* more display in case of debuging
%}
% $Date: 2013-10-31 10:05:39 +0100 (to, 31 okt 2013) $
% $Revision: 11982 $

    % Solve equations using a general iterative process defined by
    % stepFunction. This is typically Newton iterations.
    opt = struct('initialGuess',[] );
    opt = merge_options(opt, varargin{:});

    step = @(state, meta) ...
       system.stepFunction(state0, state, meta, dt, W, G, ...
                           system);

    meta = struct('converged'   , false, ...
                  'stopped'     , false, ...
                  'wellschanged', false, ...
                  'relax'       , system.nonlinear.relaxation, ...
                  'stagnate'    , false, ...
                  'iteration'   , 0);

if isempty(opt.initialGuess)
		state= state0;
   		state.wellSol = initWellSolLocal(W, state);
	else
		state = opt.initialGuess;
		state.wellSol = initWellSolLocal(W, state);
	end  

    timer = tic;
%    stagnateCount = 0;
    while ~ (meta.converged || meta.stopped),
        % Save iteration number in meta info
        meta.iteration = meta.iteration + 1;

        [state, meta,eqs] = step(state, meta);

        if meta.stagnate,
%             warning('newt:stagnate', 'Non-linear solver stagnated...')
%             stagnateCount = stagnateCount + 1;
%             if stagnateCount > 2
%                 if (~meta.converged)
%                     meta.stopped = true;
%                 end
%                     
                 break;
%             end 
%         else
%             stagnateCount = 0;
         end
    end
    if(isfield(system,'updateFinal'))
            state=system.updateFinal(state, state0);        
    end

    convergence.residuals = meta.res_history(1:meta.iteration, :);
    convergence.converged = meta.converged;
    convergence.its = meta.iteration;

    if meta.stopped,
%          warning('newt:maxit', ...
%                 ['Non-linear solver did not converge, stopped ', ...
%                  'by max iterations...']);
    end

    
    dispif(mrstVerbose, 'Completed %d iterations in %1.2f s, CNV = %e, MB = %e \n', meta.iteration, toc(timer),max(meta.CNV),max(meta.MB));
    its = meta.iteration;
end
