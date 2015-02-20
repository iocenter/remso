function [  ] = plotRunScheduleSolution( state0,solutionState,schedule )


if size(state0.s,2) == 2
    nSG = 2;
elseif size(state0.s,2) == 3
    nSG = 3;
else
    error('not supported');
end
    
    


figN = 1;

times.steps = [0;cumsum(schedule.step.val)]/day;
times.tPieceSteps = cell2mat(arrayfun(@(x)[x;x],times.steps,'UniformOutput',false));
times.tPieceSteps = times.tPieceSteps(2:end-1);


pPlot = [{state0.pressure/barsa},cellfun(@(x)x.pressure/barsa,solutionState,'UniformOutput',false)'];
sWPlot = [{state0.s(:,1)},cellfun(@(x)x.s(:,1),solutionState,'UniformOutput',false)'];
sOPlot = [{state0.s(:,2)},cellfun(@(x)x.s(:,2),solutionState,'UniformOutput',false)'];

if nSG == 3
    sGPlot = [{state0.s(:,3)},cellfun(@(x)x.s(:,3),solutionState,'UniformOutput',false)'];
    rsPlot = [{state0.rs},cellfun(@(x)x.rs,solutionState,'UniformOutput',false)'];
end



figure(figN); figN = figN+1;
plot(times.steps,cell2mat(pPlot),'-x');
ylabel('Pressure (bar)')
xlabel('time (days)')


figure(figN); figN = figN+1;
subplot(nSG,1,1);
plot(times.steps,cell2mat(sWPlot),'-x');
ylabel('Water Saturation')

subplot(nSG,1,2);
plot(times.steps,cell2mat(sOPlot),'-x');
ylabel('Oil Saturation')

if nSG == 3
    subplot(3,1,3);
    plot(times.steps,cell2mat(sGPlot),'-x');
    ylabel('Gas Saturation')
end
xlabel('time (days)')


if nSG == 3
    figure(figN); figN = figN+1;
    plot(times.steps,cell2mat(rsPlot),'-x');
    ylabel('rs (sm^3/sm^3)')
    xlabel('time (days)')
end


wellSols = cellfun(@(x)x.wellSol,solutionState,'UniformOutput',false);

if nSG ==2
    % come on MRST, there is no gas here!
    wellSols = cellfun(@(x)arrayfun(@(y)subsasgn(y,struct('type','.','subs','qGs'),0),x),wellSols,'UniformOutput',false);
end

[qWs, qOs, qGs, bhp] = wellSolToVector(wellSols);
qWs = cell2mat(arrayfun(@(x)[x,x],qWs'*day,'UniformOutput',false));
qOs = cell2mat(arrayfun(@(x)[x,x],qOs'*day,'UniformOutput',false));

if nSG == 3
    qGs = cell2mat(arrayfun(@(x)[x,x],qGs'*day,'UniformOutput',false));
end

bhp = cell2mat(arrayfun(@(x)[x,x],bhp'/barsa,'UniformOutput',false));



for ci = 1:numel(wellSols{1})
    figure(figN); figN = figN+1;
    
    subplot(nSG+1,1,1); hold all;
    plot(times.tPieceSteps, qOs(ci,:), 'x-')
    ylabel('q_o (meter^3/day)');
    title(strcat('Well: ',wellSols{1}(ci).name,' Type: ',wellSols{1}(ci).type,' Sign: ',intType2stringType(wellSols{1}(ci).sign)) );
    
    subplot(nSG+1,1,2); hold all;
    plot(times.tPieceSteps, qWs(ci,:), 'x-')
    ylabel('q_w (meter^3/day)');
    
    if nSG == 3
        subplot(4,1,3); hold all;
        plot(times.tPieceSteps, qGs(ci,:), 'x-')
        ylabel('q_G (meter^3/day)');
    end
    
    subplot(nSG+1,1,nSG+1); hold all;
    plot(times.tPieceSteps, bhp(ci,:), 'x-')
    ylabel('bhp (barsa)');
    xlabel('time(day)')
    
end



end


function stringType = intType2stringType(type)
if ischar(type)
    stringType = type;
else
    if type == 1
        stringType = 'INJ';
    elseif type == -1
        stringType = 'PROD';
    else
        warning('check converter');
        stringType = type;
    end
end

end