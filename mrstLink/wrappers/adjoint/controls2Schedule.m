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

nC = numel(schedule);
nW = numel(schedule(1).values);

uC = mat2cell(u,nW*ones(nC),1);

for i = 1:nC
   schedule(i).values = uC{i};
end


if opt.partials
    if ~isempty(opt.uScale)
        Jac = sparse(1:nu,1:nu,opt.uScale);
    else
        Jac = speye(nu);
    end
    if size(opt.uRightSeeds,1) ~= 0
        Jac = Jac*opt.uRightSeeds;
    end
    Jac = mat2cell(Jac,repmat(nW,nC,1),size(Jac,2));
else
    Jac = [];
end

end


