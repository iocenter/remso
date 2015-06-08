function [  ] = plotSolutionAdjoint(x,u,v,xd,xScale,vScale,uScale,reservoirP)

nSteps = numel(x);
qw = cell(nSteps,1);
qo = cell(nSteps,1);
bhp = cell(nSteps,1);
t0 = zeros(nSteps,1);
tf = zeros(nSteps,1);

W = reservoirP.W;
wellSys   = [ W.S ];
Dw        = blkdiag( wellSys.D );
perfSign = Dw*vertcat(W.sign);

nPerf = numel(perfSign);

if strcmp(reservoirP.S.type, 'hybrid')
   solver = 'hybrid';
else
   solver = 'mixed';
end

[state] = reservoirP.state;
for k = 1:numel(x)
    
    stepSchedule = controls2Schedule(u{reservoirP.simulationControlSteps(k)},reservoirP.schedule(k),'uScale',uScale);
    
    W     = updateWells(reservoirP.W, stepSchedule);
    
    resSol = solveIncompFlowLocal(state, reservoirP.G, reservoirP.S, reservoirP.fluid, ...
        'gravityOff',false,...
        'wells', W, 'Solver', solver);
    
    simRes.resSol = resSol;
    simRes.wellSol = resSol.wellSol;
    simRes.timeInterval = stepSchedule.timeInterval;
    
	[state] = stateVector2stateMrst( x{k},reservoirP,'xScale',xScale);

    simRes.resSol.s = state.s;
    
    [qw(k),qo(k),bhp(k),t0(k),tf(k)] = perfVariables(W, reservoirP.fluid, simRes);
    
    if norm((qw{k}+qo{k}).*(perfSign)-v{k}(1:nPerf).*vScale(1:nPerf)) > sqrt(eps)
        warning('The solution seems not correct');
    end

end


t = [t0,tf];
t = reshape(t',numel(t),1);
qwP = cell2mat(cellfun(@(q)[q,q]',qw,'UniformOutput',false));
qoP = cell2mat(cellfun(@(q)[q,q]',qo,'UniformOutput',false));
bhpP = cell2mat(cellfun(@(q)[q,q]',bhp,'UniformOutput',false));

figure(1);plot(t/day,qwP*day)
figure(2);plot(t/day,qoP*day)
figure(3);plot(t/day,bhpP/barsa)

end

