function [schedule] = eclipseSchedule2mrstSchedule(schedule,G,rock)
%
% Transform a eclipse schedule to a mrst schedule
%
useMrstSchedule = isfield(schedule.control(1), 'W');

if ~useMrstSchedule
    nc = numel(schedule.control);
    mrstControl = repmat(struct('W',[]),nc,1);
    for k = 1:nc
        % processWells is call as in runScheduleADI
        mrstControl(k).W =  processWells(G, rock, schedule.control(k),...
            'Verbose', true, 'DepthReorder', false);
    end
    schedule.control = mrstControl;
    
end


end