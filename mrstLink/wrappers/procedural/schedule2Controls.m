function [ u,Jac] = schedule2Controls(schedule,varargin)
%
% extract the schedule control information and scale it
%
%

opt = struct('uScale',[]);
opt = merge_options(opt, varargin{:});

Jac = [];

[ vals ] = schedule2CellControls(schedule);
u = cellControls2Controls(vals);

if ~isempty(opt.uScale)
    u = u./opt.uScale;
    
    if nargout > 1
        nu = numel(u);
        Jac = sparse(1:nu,1:nu,1./opt.uScale);
    end
else
    if nargout > 1
        nu = numel(u);
        Jac = speye(nu);
    end
end


end