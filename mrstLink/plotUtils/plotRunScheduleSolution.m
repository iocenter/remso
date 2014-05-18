function [  ] = plotRunScheduleSolution( solutionState,schedule )


figN = 1;

times.steps = [0;cumsum(schedule.step.val)]/day;
times.tPieceSteps = cell2mat(arrayfun(@(x)[x;x],times.steps,'UniformOutput',false));
times.tPieceSteps = times.tPieceSteps(2:end-1);


pPlot = cellfun(@(x)x.pressure/barsa,solutionState,'UniformOutput',false)';
sPlot = cellfun(@(x)x.s(:,1),solutionState,'UniformOutput',false)';



figure(figN); figN = figN+1;
plot(times.steps,cell2mat(pPlot),'-x');
ylabel('Pressure (bar)')
xlabel('time (days)')


figure(figN); figN = figN+1;
plot(times.steps,cell2mat(sPlot),'-x');
ylabel('Saturation')
xlabel('time (days)')



wellSols = cellfun(@(x)x.wellSol,solutionState,'UniformOutput',false);
wellSols = wellSols(2:end);

[qWs, qOs, qGs, bhp] = wellSolToVector(wellSols);
qWs = cell2mat(arrayfun(@(x)[x,x],qWs'*day,'UniformOutput',false));
qOs = cell2mat(arrayfun(@(x)[x,x],qOs'*day,'UniformOutput',false));
bhp = cell2mat(arrayfun(@(x)[x,x],bhp'/barsa,'UniformOutput',false));


figure(figN); figN = figN+1;
for ci = 1:size(wellSols{1},2)
    
    
    subplot(3,1,1); hold all;
    plot(times.tPieceSteps, qOs(ci,:), 'x-')
    ylabel('q_o (meter^3/day)');
    title(strcat('Well: ',wellSols{1}(ci).name,' Type: ',wellSols{1}(ci).type,' Sign: ',intType2stringType(wellSols{1}(ci).sign)) );
    
    subplot(3,1,2); hold all;
    plot(times.tPieceSteps, qWs(ci,:), 'x-')
    ylabel('q_w (meter^3/day)');
    
    
    subplot(3,1,3); hold all;
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