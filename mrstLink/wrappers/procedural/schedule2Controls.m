function [ u,Jac] = schedule2Controls(schedule,varargin)
%
% extract the schedule control information and scale it
%
%

opt = struct('uScale',[], 'fixedWells', []);
opt = merge_options(opt, varargin{:});

Jac = [];

[ valsFull ] = schedule2CellControls(schedule);


nW = numel(schedule.control(1).W);
controlWells = setdiff(1:nW, opt.fixedWells);
assert(all(arrayfun(@(c) numel(c.W),  schedule.control)==nW)); %% all controls have the same number of wells

vals = cellfun(@(v) v(controlWells), valsFull, 'UniformOutput', false);  %% remove fixed wells

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