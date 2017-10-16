


schedule = reservoirP.schedule;

simulationSteps = repmat(5,365/5*10,1); 

controlSteps   =    [82,24,24,24,12,12,6 ,4  ,2  ,1];
controlStepLength = [5 ,15,15,15,30,30,60,90,180,360];
controlStepRep =    controlStepLength/5;
repmat(5,365/5*10,1)


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
 
         
