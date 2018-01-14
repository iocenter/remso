function ok = testTargetGradient(state,model,schedule,target,varargin)

if numel(varargin)>=1
    reducedTest=varargin{1};    
else
    reducedTest=true;
end
if numel(varargin)>=2
    pert=varargin{2};    
else
    pert=1e-5;
end
if numel(varargin)>=3
    tol=varargin{3};    
else
    tol=1e-3;
end


oks = [];

[wellSols, states, report] = simulateScheduleAD(state, model, schedule);

uScale = schedule2Controls(schedulesScaling(schedule,...
    'RATE',model.scaling.qWs,...
    'ORAT',model.scaling.qOs,...
    'WRAT',model.scaling.qWs,...
    'LRAT',model.scaling.qWs,...
    'RESV',0,...
    'BHP',model.scaling.bhp));


% test transformations
[u, J1] = schedule2Controls(schedule,'uScale',uScale);
ut = rand(size(u));
[ schedulet, J2 ] = controls2Schedule(ut,schedule,'uScale',uScale,'partials',true);
[u2,~] = schedule2Controls(schedulet,'uScale',uScale);
uDims = cellfun(@(ui)numel(ui),schedule2CellControls(schedule));

uErr = norm(ut-u2,inf);
ugErr = norm(J1*vertcat(J2{:})-eye(numel(u)),inf);
oks = horzcat(oks,uErr< tol); %#ok<*AGROW>
oks = horzcat(oks,ugErr< tol);


nx = model.G.cells.num*(sum(model.getActivePhases()));

seeds = speye(nx+sum(uDims));
xRightSeeds = seeds(1:nx,:);
uRightSeeds = mat2cell(seeds(nx+1:end,:),uDims,nx+sum(uDims));

objective = @(tstep) target(wellSols, schedule, 'ComputePartials', true, 'tStep', tstep);

o1 = objective(1);
nt = numel(double(o1{1}));

forwardGradient = computeGradientForwardAD(state, states, model,schedule, objective,xRightSeeds,uRightSeeds);
adjointGradient = computeGradientAdjointAD(state, states, model,schedule, objective,'xRightSeeds',xRightSeeds,'uRightSeeds',uRightSeeds);

absError = @(a,b) norm(a-b,inf);

errAdjFwd = absError(forwardGradient,adjointGradient);

adjGradX0 = adjointGradient(:,1:nx);
adjGradU = adjointGradient(:,nx+1:end);


[x0,dXdState] = model.toStateVector(state);
[u,dUdSchedule]  = schedule2Controls(schedule,'uScale',uScale);


[~,dStatedX] = model.toMRSTStates(x0);
[~,dScheduledU]  = controls2Schedule(u,schedule,'uScale',uScale,'partials',true);

if reducedTest
    pu = rand(numel(u),1).*abs(u);
else
    pu = diag(abs(u));
end
grad = zeros(nt,size(pu,2));
for k = 1:size(pu,2)
    uP = u + pu(:,k)*pert;
    uN = u - pu(:,k)*pert;    

    grad(:,k) = (testObjective(target,model,schedule,x0,uP,uScale ) - testObjective(target,model,schedule,x0,uN,uScale ))/pert/2;
end

scaledAdjGrad = (adjGradU*vertcat(dScheduledU{:}));
errU = absError(grad,scaledAdjGrad*pu);

if reducedTest
    d = rand(numel(x0),1);
    px = d.*abs(x0)/norm(d);
else
    px = diag(abs(x0));
end    
gradx0 = zeros(nt,size(px,2));
for k = 1:size(px,2)
    x0P = x0 + px(:,k)*pert;
    x0N = x0 - px(:,k)*pert;

    gradx0(:,k) = (testObjective(target,model,schedule,x0P,u,uScale ) - testObjective(target,model,schedule,x0N,u,uScale ))/pert/2;
end

errX = absError(gradx0,(adjGradX0*dStatedX*px));

fprintf('Maximum absolute error in control conversion = %e\n', uErr)
fprintf('Maximum absolute error in control conversion jacobian = %e\n', ugErr)
fprintf('Maximum absolute gradient error between AD methods = %e\n',errAdjFwd)
fprintf('Maximum absolute gradient error w.r.t controls = %e\n',errU)
fprintf('Maximum absolute gradient error w.r.t states = %e\n',errX)

oks = horzcat(oks,errAdjFwd<tol,errU<tol,errX<tol);

ok = all(oks);

end

function [ objective ] = testObjective(target,model,schedule,state,u,uScale)

[ state] = model.toMRSTStates( state);
[ schedule ] = controls2Schedule( u,schedule,'uScale',uScale);

[wellSols, states, report] = simulateScheduleAD(state, model, schedule);

objective = target(wellSols, schedule);

objective = sum(cat(3,objective{:}),3);


end
