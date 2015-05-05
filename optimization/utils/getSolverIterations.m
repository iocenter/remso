function [its] = getSolverIterations(simVars)
% get the number of iterations!


if isa(simVars,'Composite')
    spmd
        itsC = cellfun(@getConv,simVars);
    end


    its = cell(numel(itsC),1);
    for k = 1:numel(itsC)
        its{k} = itsC{k};
    end
    its = cell2mat(its);
else
    its = cellfun(@getConv,simVars);
end



end

function c = getConv(sV)
    if iscell(sV)
        c = mean(cellfun(@getConv,sV));
    else
        c = sV.convergence.its;
    end
end
