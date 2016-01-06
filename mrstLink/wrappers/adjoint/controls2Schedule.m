function [ schedule,Jac ] = controls2Schedule( u,schedule,controls,W,varargin)
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
nW = numel(controls(1).well);

uC = mat2cell(u,nW*ones(nC),1);


[A_N, b_N, A_D, b_D] = controls2Wells(W, schedule, controls);


rateWells     =  strcmp('rate', {W.type}) ;
BHPWells      =  strcmp('bhp', {W.type}) ;

for n = 1:nC
    
	q = A_N{n}*uC{n} + b_N{n};
	p = A_D{n}*uC{n} + b_D{n};    
    
	schedule(n).values(rateWells) = q;
    schedule(n).values(BHPWells) = p;

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


