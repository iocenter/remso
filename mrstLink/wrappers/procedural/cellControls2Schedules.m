function [ schedules,Jac] = cellControls2Schedules(uCells,schedules,varargin)
%
%  write the cellcontrols in the schedule form
%
%
opt = struct('cellControlScales',[]);
opt = merge_options(opt, varargin{:});


partials = (nargout>1);
n_sp = numel(schedules);
Jac= cell(n_sp,1);


for k = 1:n_sp
    if ~isempty(opt.cellControlScales)
        [schedules(k),Jac{k}] = controls2Schedule(uCells{k},schedules(k),'uScale',opt.cellControlScales{k},'partials',partials);
    else
        [schedules(k),Jac{k}] = controls2Schedule(uCells{k},schedules(k),'partials',partials);
    end
end




end

