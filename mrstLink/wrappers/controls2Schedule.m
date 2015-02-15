function [ schedule,Jac ] = controls2Schedule( u,schedule,varargin)
%
%  set the controls in the schedule.  Scale the variables
%
%

opt = struct('uScale',[],'partials',false,'uRightSeeds',[]);
opt = merge_options(opt, varargin{:});

nu = numel(u);

if ~isempty(opt.uScale)
    u = u.*opt.uScale;
end

[vals,nC,nW] = controls2CellControls(u,schedule);

schedule = cellControls2Schedule(vals,schedule);


if opt.partials
    if ~isempty(opt.uScale)
        Jac = bsxfun(@times,speye(nu),opt.uScale');
    else
        Jac = speye(nu);
    end
    if size(opt.uRightSeeds,1) ~= 0
        Jac = Jac*opt.uRightSeeds;
    end
    Jac = mat2cell(Jac,nW*ones(nC,1),size(Jac,2));
else
    Jac = [];
end

end


