function [ target,lbt,ubt,sparsity] = stochasticOutputVarsSelector(lbs,ubs,uDims)
%
% Create a function that returns the state constraints x and v with finite
% bounds
%



finiteBoundsS = or(isfinite(lbs),isfinite(ubs));

lbt = lbs(finiteBoundsS);
ubt = ubs(finiteBoundsS);


sparsity = ones(size(lbt,1),sum(uDims)); %% can we do better without knowing what is S?
% this shouldn't be too big anyway


target = @(s,u,varargin)outputSelectionS(s,u,finiteBoundsS,varargin{:});



end