function [vals] = controls2CellControls(u,schedule)
%
%  divide the controls in cells according to the number of wells
%

useMrstSchedule = isfield(schedule.control(1), 'W');

% assuming same number of wells during the schedule
if useMrstSchedule
    numC = numel(schedule.control);
    numW  = numel(schedule.control(1).W);
else
    numC = numel(schedule.control);
    numW  = size(schedule.control(1).WCONINJE, 1)+size(schedule.control(1).WCONPROD, 1);
end

vals = mat2cell(u,numW*ones(numC,1),1);



end