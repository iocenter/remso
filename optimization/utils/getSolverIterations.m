function [its] = getSolverIterations(simVars)
% get the number of iterations!


if isa(simVars,'Composite')
    spmd
        itsC = cellfun(@getConv,simVars);
    end


    its = [];
    for k = 1:numel(itsC)
        its = [its;itsC{k}];
    end
else
    its = cellfun(@getConv,simVars);
end



end

function c = getConv(sV)
    c = sV.convergence.its;
end
