function [ schedule ] = cellControls2Schedule( vals,schedule )

%
% write the cell controls in the schedule
%

useMrstSchedule = isfield(schedule.control(1), 'W');

if useMrstSchedule
    
    for k = 1:numel(schedule.control)
        [ schedule.control(k).W ] = values2Wells( vals{k},schedule.control(k).W );
        
    end
    
else
    schedule = updateSchedule(schedule, vals);  %% cellcontrols2Schedule
end

%--------------------------------------------------------------------------
end



