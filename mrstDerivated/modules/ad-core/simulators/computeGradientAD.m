function grad = computeGradientAD(state0, states, model, schedule, objective,xRightSeeds,uRightSeeds,varargin)
% From https://github.com/iocenter/remso/ 4f8fa54a92c14b117445ad54e9e5dd3a0e47a7f5

%TODO: how to decide on the linear solver?

nAdj = size(objective(1),1);
nFwd = size(xRightSeeds,2);


if (nFwd > nAdj )
    grad = computeGradientAdjointAD(state0, states, model, schedule, objective,'xRightSeeds',xRightSeeds,'uRightSeeds',uRightSeeds,'LinearSolver',[]);
else
    grad = computeGradientForwardAD(state0, states, model, schedule, objective,xRightSeeds,uRightSeeds,'LinearSolver',[]);
end


end