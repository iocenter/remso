function [ schedules ] = cellControls2Schedules(uCells,schedules,varargin)
%
%  write the cellcontrols in the schedule form
%
%
opt = struct('cellControlScales',[]);
opt = merge_options(opt, varargin{:});



n_sp = numel(schedules);

for k = 1:n_sp
    if ~isempty(opt.cellControlScales)
        schedules(k) = controls2Schedule(uCells{k},schedules(k),'uScale',opt.cellControlScales{k});
    else
        schedules(k) = controls2Schedule(uCells{k},schedules(k));
    end
end




end

