function [ stateMrst,Jac ] = stateVector2stateMrst( stateVector,reservoirP,varargin)
%
%  write a state vector as a mrst state structure
%
%

opt = struct('xScale',[],'partials',false);
opt = merge_options(opt, varargin{:});


if ~isempty(opt.xScale)
    stateVector = stateVector.*opt.xScale;
end



stateMrst = initResSol(reservoirP.G, 0.0,stateVector);
stateMrst.wellSol = initWellSol(reservoirP.W, 0);

if opt.partials
    nx = numel(stateVector);

    if ~isempty(opt.xScale)
        Jac = bsxfun(@times,speye(nx),opt.xScale');
    else
        Jac = speye(nx);
    end
else
    Jac = [];
end




end