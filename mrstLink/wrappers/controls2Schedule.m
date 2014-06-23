function [ schedule,Jac ] = controls2Schedule( u,schedule,varargin)
%
%  set the controls in the schedule.  Scale the variables
%
%

opt = struct('uScale',[],'partials',false);
opt = merge_options(opt, varargin{:});

nu = numel(u);

if ~isempty(opt.uScale)
    u = u.*opt.uScale;
end

vals = controls2CellControls(u,schedule);

schedule = cellControls2Schedule(vals,schedule);


if opt.partials
    if ~isempty(opt.uScale)
        Jac = bsxfun(@times,speye(nu),opt.uScale');
    else
        Jac = speye(nu);
    end
else
    Jac = [];
end

end


