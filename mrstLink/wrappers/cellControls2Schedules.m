function [ schedules ] = cellControls2Schedules(uCells,schedules,varargin)
%
%  write the cellcontrols in the schedule form
%
%
opt = struct('doScale',false,'cellControlScales',[]);
opt = merge_options(opt, varargin{:});



n_sp = numel(schedules);

for k = 1:n_sp
    if opt.doScale
        schedules(k) = controls2Schedule(uCells{k},schedules(k),'doScale',true,'uScale',opt.cellControlScales{k});
    else
        schedules(k) = controls2Schedule(uCells{k},schedules(k),'doScale',false);
    end
end




end

