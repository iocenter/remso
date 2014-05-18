function [ violation ] = checkConstraintViolation(x,lb,lowActive,ub,upActive )
%  violation = sum of all the constraints violation (norm1 of the violation)


% if you want infinity norm
% geqZeroNorm = @(x)max([x;0]);
%  violation = max([cellfun(@(var,upBound,act) geqZeroNorm(var(act) -upBound(act)),x,ub,upActive);...
%                   cellfun(@(var,lowBound,act)geqZeroNorm(lowBound(act)-var(act)),x,lb,lowActive)]);


geqZeroNorm = @(x)sum(max([x;0],0));

violation = sum([cellfun(@(var,upBound,act) geqZeroNorm(var(act) -upBound(act)),x,ub,upActive);...
                 cellfun(@(var,lowBound,act)geqZeroNorm(lowBound(act)-var(act)),x,lb,lowActive)]);



end

