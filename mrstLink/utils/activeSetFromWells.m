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
activeStatesV = wellSol2algVar(activeStatesV);
activeStatesV = true(size(activeStatesV));



injIndex = vertcat(W.sign) <0;

wInj = W(injIndex);
wProd = W(~injIndex);


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

TODO: create a function to instantiate activeSet

[lowActiveS] = stateMrst2stateVector( lowActiveS );
lowActive.x = repmat({lowActiveS},totalPredictionSteps,1);
lowActive.v = repmat({[activeStatesV;true]},totalPredictionSteps,1);

[upActiveS] = stateMrst2stateVector( upActiveS );
upActive.x = repmat({upActiveS},totalPredictionSteps,1);
upActive.v = repmat({[activeStatesV;true]},totalPredictionSteps,1);

