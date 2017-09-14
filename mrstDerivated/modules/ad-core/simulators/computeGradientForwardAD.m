function gradients = computeGradientForwardAD(state0, states, model, schedule, objective,xRightSeeds,uRightSeeds,varargin)
% From https://github.com/iocenter/remso/ 4f8fa54a92c14b117445ad54e9e5dd3a0e47a7f5
%Compute gradients using an forward simulation that is linear in each step


    opt = struct('LinearSolver',[]);
    opt = merge_options(opt, varargin{:});
    
    getState = @(i) getStateFromInput(schedule, states, state0, i);
    
    if isempty(opt.LinearSolver)
        linsolve = BackslashSolverAD();
    else
        linsolve = opt.LinearSolver;
    end
      
    nstep = numel(schedule.step.val);
    nt = nstep;
    
    gradients = 0;
    grad = xRightSeeds;
    for step = 1:nt
        uRightSeed = uRightSeeds{schedule.step.control(step)};
        [grad,report] = model.solveForwardGrad( linsolve, getState,...
                 schedule,step,grad,uRightSeed);
        gradients = gradients + extractJac(objective(step))*grad;
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

function jac = extractJac(objective)
if ~isnumeric(objective)
    if iscell(objective)
        objective = objective{:};
    end
    objective = cat(objective);
    
    % Above CAT means '.jac' is a single element cell array.  Extract contents.
    jac  = objective.jac{1};
else
    % the actual jacobian is given
    jac  = objective;
end
end
