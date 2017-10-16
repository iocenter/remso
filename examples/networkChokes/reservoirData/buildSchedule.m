d = 5;


schedule = reservoirP.schedule;

simulationSteps = repmat(d,(365/d)*10,1); 

controlSteps   =    [82,24,24,24,12,12,6 ,4  ,2  ,1];

% controlSteps = [24, 12, 12, 6, 6, 3, 3, 2, 1, 1];

controlStepLength = [d ,3*d,3*d,3*d,6*d,6*d,12*d,18*d,36*d,72*d];
controlStepRep =    controlStepLength/d;
repmat(d,(365/d)*10,1);


controlIndex = [];
ic = 1;
for i = 1:numel(controlSteps)
    for k = 1:controlSteps(i)
        controlIndex = [controlIndex;repmat(ic,controlStepRep(i),1)];
        ic = ic+1;
    end
    
end

reservoirP.schedule.step.control = controlIndex;
reservoirP.schedule.step.val = simulationSteps*day;
reservoirP.schedule.control = repmat(reservoirP.schedule.control(1),max(controlIndex),1);
 
         
