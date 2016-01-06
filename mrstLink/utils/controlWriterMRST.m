function [] = controlWriterMRST(u,it,controlSchedules,cellControlScales,varargin)

opt = struct('fileName', 'schedule.inc','units','METRIC');
opt = merge_options(opt, varargin{:});


schedules = cellControls2Schedules(u,controlSchedules,'cellControlScales',cellControlScales);
writeSchedule(schedules,'fileName',opt.fileName,'units',opt.units);


end

