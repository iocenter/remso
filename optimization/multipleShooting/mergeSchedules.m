function schedule = mergeSchedules(schedules)
%
%  Join the separated schedules into a single schedule
%
%

n_sp = numel(schedules);

schedule = schedules(1);

for k = 2:n_sp
    
    nc = numel(schedule.control);
    
    schedule.control = [schedule.control;schedules(k).control]; 
    schedule.step.control = [schedule.step.control;nc+schedules(k).step.control];
    schedule.step.val = [schedule.step.val;schedules(k).step.val];
    
end
schedule.time = schedules(1).time;


end