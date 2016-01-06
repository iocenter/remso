function grad = computeGradientAD(state0, states, model, schedule, objective,xRightSeeds,uRightSeeds,varargin)

%TODO: how to decide on the linear solver?

nAdj = size(objective(1),1);
nFwd = size(xRightSeeds,2);


if (nFwd > nAdj )
    grad = computeGradientAdjointAD(state0, states, model, schedule, objective,'xRightSeeds',xRightSeeds,'uRightSeeds',uRightSeeds,'LinearSolver',[]);
else
    grad = computeGradientForwardAD(state0, states, model, schedule, objective,xRightSeeds,uRightSeeds,'LinearSolver',[]);
end


end