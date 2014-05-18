function [ nW ] = getnW( schedule )
%
% given a schedule, provide the number of wells

useMrstSchedule = isfield(schedule.control(1), 'W');
if useMrstSchedule
    nW = numel(schedule.control(1).W);
else
    nW = size(schedule.control(1).WELSPECS,1);    
end



end

