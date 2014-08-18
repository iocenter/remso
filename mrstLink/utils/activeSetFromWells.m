function [ lowActive,upActive ] = activeSetFromWells( reservoirP,totalPredictionSteps)
%
%  Generate a guess of the active sets.  Consider all well variables and
%  grid-blocks with perforation variables in the first QP



if (isfield(reservoirP.schedule.control,'W'))
   W =  reservoirP.schedule.control.W;
else
   W = processWellsLocal(reservoirP.G, reservoirP.rock,reservoirP.schedule.control(1),'DepthReorder', true);
end


activeStatesV = initWellSolLocal(W, reservoirP.state);
activeStatesV = wellSol2algVar(activeStatesV,...
    'activeComponents',reservoirP.system.activeComponents);
activeStatesV = true(size(activeStatesV));



injIndex = vertcat(W.sign) <0;

wInj = W(injIndex);
wProd = W(~injIndex);



comp = reservoirP.system.activeComponents;
if ~comp.gas && ~comp.polymer && ~(comp.T || comp.MI)
    
    activeStates = reservoirP.state;
    activeStates.pressure = false(size(activeStates.pressure));
    activeStates.s = false(size(activeStates.s));
    
    lowActiveS = activeStates;
    upActiveS = activeStates;
    
    wcInj = vertcat(wInj.cells);
    for k = wcInj
        upActiveS.pressure(k) = true;
        lowActiveS.s(k,1) = true(size(activeStates.s(k,1)));
        upActiveS.s(k,1) = true(size(activeStates.s(k,1)));
        lowActiveS.pressure(k) = true;
    end
    
    wcProd = vertcat(wProd.cells);
    for k = wcProd
        lowActiveS.pressure(k) = true;
        upActiveS.s(k,1) = true(size(activeStates.s(k,1)));
        upActiveS.pressure(k) = true;
        lowActiveS.s(k,1) = true(size(activeStates.s(k,1)));
    end
        
    
    
    [lowActiveS] = stateMrst2stateVector( lowActiveS);
    lowActive.x = repmat({lowActiveS},totalPredictionSteps,1);
    
    [upActiveS] = stateMrst2stateVector( upActiveS);
    upActive.x = repmat({upActiveS},totalPredictionSteps,1);
    
    
    
elseif comp.gas && comp.oil && comp.water
    
    
    activeState = false(reservoirP.G.cells.num,1);
    cells = [vertcat(wInj.cells);vertcat(wProd.cells)];
    for k = cells
        activeState(k) = true;
    end
    activeStates = repmat(activeState,3,1);

    lowActive.x = repmat({activeStates},totalPredictionSteps,1);
    upActive.x = repmat({activeStates},totalPredictionSteps,1);

else
    error('Not implemented for current activeComponents');
end


lowActive.v = repmat({activeStatesV},totalPredictionSteps,1);
upActive.v = repmat({activeStatesV},totalPredictionSteps,1);

