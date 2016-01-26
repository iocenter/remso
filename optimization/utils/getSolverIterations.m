function [its] = getSolverIterations(simVars)
% get the number of iterations!


if isa(simVars,'Composite')
    spmd
        itsC = cell2mat(cellfun(@getConv,simVars,'UniformOutput',false));
    end


    its = zeros(numel(itsC),1);
    for k = 1:numel(itsC)
        its(k) = itsC{k};
    end
else
    its = cell2mat(cellfun(@getConv,simVars,'UniformOutput',false));
end



end

function c = getConv(sV)
    if iscell(sV)
        c = cellfun(@getConv,sV);
    else
        c = sV.convergence.its;
    end
end
