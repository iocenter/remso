function [state, dx, fail,i] = linesearchADI(state0, dx0, system, getEqs, updateState, isBO)
%{
Modification by codas

% use norm1 only.
% change backtrack parameter to 5
% modify acceptance condition to armijo rule

%}
    ni = system.nonlinear.lineIter;
    target = system.nonlinear.lineRelTol;

    
    fail = true;
    i = 0;
    alph = 0;

    e = @(eqs) cellfun(@(x) norm(x, 1), {eqs{system.cellwise}});

    if isBO
        [eqs, history, explTrms] = getEqs(state0);
    else
        eqs = getEqs(state0);
    end
    err0 = e(eqs);


    dx = dx0;

    % obs!
    % error0 is the error at 0 and -error0 its gradient on the line dx
    error0 = norm(err0,1);
    lineTarget = (1-target)*error0;
    while fail && (i < ni)
        dx = cellfun(@(x) x* 5^alph, dx0, 'UniformOutput', false);
        if isBO
            state = updateState(dx, explTrms);
            [eqs, history, explTrms] = getEqs(state);
        else
            state = updateState(dx);
            eqs = getEqs(state);
        end
        err = e(eqs);

        
        fail = norm(err,1) > error0 - lineTarget*5^(alph);

        alph = alph - 1;
        i    = i + 1;

    end
    if fail
        state = state0;
        dx = dx0;
    end
end
