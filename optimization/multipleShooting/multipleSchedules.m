function schedules = multipleSchedules(schedule, shootingIntervals)
%
% Divide an MRST schedule in several schedules according to the shooting
% Intervals
%
%

n_sp = numel(shootingIntervals);

schedules = repmat(schedule,n_sp, 1);

times0 = cumsum(schedule.step.val) + schedule.time;


indexP = 1;

for k = 1:n_sp
    
    step = struct('control',schedule.step.control(indexP:shootingIntervals(k)),'val',schedule.step.val(indexP:shootingIntervals(k)));
    
    schedules(k).control = schedule.control;
    schedules(k).step =step;
    if(k==1)        
        schedules(1).time = schedule.time;
    else
        schedules(k).time = times0(shootingIntervals(k-1));
    end
    
    schedules(k) = removeUnreferencedControls(schedules(k));
    
    indexP = shootingIntervals(k)+1;
end

end
function schedule = removeUnreferencedControls(schedule)

n = numel(schedule.step.control);

found = zeros(1,n);

controlValues = schedule.control;  % memory allocating
controlReference = zeros(n,1);

c = 0;
for i = 1:n
    
    if(found(i))
        continue;
    end
    c = c + 1;
    
    f = find(schedule.step.control==schedule.step.control(i));
    found(f) = found(f) + 1;
    
    controlValues(c) = schedule.control(schedule.step.control(i));
    
    for j = f
        controlReference(j) = c;
    end
end

schedule.control = controlValues(1:c);
schedule.step.control = controlReference;

assert(prod(found) == 1);

end