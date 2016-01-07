function [] = controlWriterMRST(u,it,controlSchedules,cellControlScales,varargin)

opt = struct('fileName', 'schedule.inc','units','METRIC', 'fixedWells', []);
opt = merge_options(opt, varargin{:});


schedules = cellControls2Schedules(u,controlSchedules,'cellControlScales',cellControlScales, 'fixedWells', opt.fixedWells);
writeSchedule(schedules,'fileName',opt.fileName,'units',opt.units);


end

