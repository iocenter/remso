function [ schedule ] = controls2Schedule( u,schedule,varargin)
%
%  set the controls in the schedule.  Scale the variables
%
%

opt = struct('doScale',false,'uScale',[]);
opt = merge_options(opt, varargin{:});



if opt.doScale
    if isempty(opt.uScale)
        [ scheduleScale ] = schedulesScaling( schedule );
        [ uScale ] = schedule2Controls(scheduleScale);
        u = u.*uScale;
    else
        u = u.*opt.uScale;
    end
end

vals = controls2CellControls(u,schedule);

schedule = cellControls2Schedule(vals,schedule);  

