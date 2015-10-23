function [  ] = plotSolutionOW( x,u,v,d, lbv, ubv, lbu, ubu, ss,obj,times,xScale,uScale,vScale, nScale, uScalePlot,schedules,wellSol,lbuPot,ubuPlot,ulim,minState,maxState,varargin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% varargin = {'simulate',[],'xScale',xScale,'uScale',cellControlScales,'uScalePlot',cellControlScalesPlot,'schedules',mShootingP.schedules}
% opt = struct('simulate',[],'simFlag',false,'plotWellSols',true,'plotSchedules',true,'plotObjective',true,'pF',@(x)x,'sF',@(x)x,'figN',1000,'wc',false,'reservoirP',[],'plotSweep',false);
opt = struct('simulate',[],'simFlag',false,'plotWellSols',true, 'plotNetsol', true, 'plotNetControls', true, 'plotSchedules',true,'plotObjective',true,'pF',@(x)x,'sF',@(x)x,'figN',1000,'wc',false,'reservoirP',[],'plotSweep',false);
opt = merge_options(opt, varargin{:});

figN = opt.figN;

xM = cellfun(@(xi)stateVector2stateMrst( xi,'xScale',xScale),[ss.state;x],'UniformOutput',false);
dM = cellfun(@(xi)stateVector2stateMrst( xi,'xScale',xScale),d,'UniformOutput',false);


pPlot = cellfun(@(x)opt.pF(x.pressure/barsa),xM,'UniformOutput',false)';
sPlot = cellfun(@(x)opt.sF(x.s(:,1)),xM,'UniformOutput',false)';
dpPlot = cellfun(@(x)opt.pF(x.pressure/barsa),dM,'UniformOutput',false)';
dsPlot = cellfun(@(x)opt.sF(x.s(:,1)),dM,'UniformOutput',false)';


figure(figN); figN = figN+1;
subplot(2,1,1)
plot(times.steps,cell2mat(pPlot),'-x');
ylabel('Pressure (bar)')
xlabel('time (days)')
ylim([minState.pressure,maxState.pressure]/barsa)
subplot(2,1,2)
plot(times.steps(2:end),cell2mat(dpPlot),'-x');
ylabel('Pressure error (bar)')
xlabel('time (days)')



figure(figN); figN = figN+1;
subplot(2,1,1)
plot(times.steps,cell2mat(sPlot),'-x');
ylabel('Saturation')
xlabel('time (days)')
ylim([minState.sW,maxState.sW])
subplot(2,1,2)
plot(times.steps(2:end),cell2mat(dsPlot),'-x');
ylabel('Saturation error')
xlabel('time (days)')


if ~isempty(opt.reservoirP) && opt.plotSweep
    figure(figN); figN = figN+1;
    for k = 1:numel(x)
        plotCellData(opt.reservoirP.G,xM{k}.s(:,1));
        title(['Saturation - year: ' num2str(times.steps(k)/year)])
        colorbar
        pause(0.01)
    end
end

nW = arrayfun(@(s)numel(s.control.W),schedules,'UniformOutput',false);

w = cellfun(@(ui,nw)ui(1:nw),u,nW,'UniformOutput',false);
p = cellfun(@(ui,nw)ui(nw+1:end),u,nW,'UniformOutput',false);

wScale = cellfun(@(ui,nw)ui(1:nw),uScale,nW,'UniformOutput',false);
wScalePlot = cellfun(@(ui,nw)ui(1:nw),uScalePlot,nW,'UniformOutput',false);

pScale = cellfun(@(ui,nw)ui(nw+1:end),uScale,nW,'UniformOutput',false);

if opt.plotSchedules
    [uM,schedulesSI] = scaleSchedulePlot(w,schedules,wScale,wScalePlot);
    
    uPiece = cell2mat(arrayfun(@(x)[x,x],uM,'UniformOutput',false));
    
    [~,type,control] =schedule2CellControls(schedules(1));
    
    type = cellfun(@(x)intType2stringType(x),type{1},'UniformOutput',false);
    
    for ci = 1:size(uPiece,1)
        figure(figN); figN = figN+1;
        plot(times.tPieceControls,uPiece(ci,:),'bo-',times.tPieceControls,lbuPot(ci,:),'rx-',times.tPieceControls,ubuPlot(ci,:),'rx-')
        ylim([ulim(ci,1) ulim(ci,2)])
        legend([num2str(ci),' ',type{ci},' ',control{1}{ci}])
        xlabel('time (days)')
        
    end
end


if opt.plotWellSols
    [uM,schedulesSI] = scaleSchedulePlot(w,schedules,wScale,wScalePlot);
    
    
    totalPredictionSteps = getTotalPredictionSteps(ss);
    wellSols = cell(1,totalPredictionSteps);
    for k = 1:totalPredictionSteps
        wellSols{k} = algVar2wellSol(v{k},wellSol,'vScale',vScale);
    end
    
    
    %plot results
    
    [qWs, qOs, qGs, bhp] = wellSolToVector(wellSols);
    qWs = cell2mat(arrayfun(@(x)[x,x],qWs'*day,'UniformOutput',false));
    qOs = cell2mat(arrayfun(@(x)[x,x],qOs'*day,'UniformOutput',false));
    bhp = cell2mat(arrayfun(@(x)[x,x],bhp'/barsa,'UniformOutput',false));
    
    simulateFlag = ~isempty(opt.simulate) && opt.simFlag;
    if simulateFlag
        
        
        initialGuess = cell(1,totalPredictionSteps);
        for k = 1:totalPredictionSteps
            initialGuess{k} = xM{k};
            initialGuess{k}.wellSol = wellSols{k};
        end
        
        schedule = mergeSchedules(schedulesSI);
        [wellSolsS,~,scheduleSol] = opt.simulate(schedule,'initialGuess',initialGuess);
        [qWsS, qOsS, qGsS, bhpS] = wellSolToVector(wellSolsS);
        qWsS = cell2mat(arrayfun(@(x)[x,x],qWsS'*day,'UniformOutput',false));
        qOsS = cell2mat(arrayfun(@(x)[x,x],qOsS'*day,'UniformOutput',false));
        bhpS = cell2mat(arrayfun(@(x)[x,x],bhpS'/barsa,'UniformOutput',false));
        
        timesSol.steps = (cumsum([scheduleSol.time;scheduleSol.time + scheduleSol.step.val]))/day;
        timesSol.tPieceSteps = cell2mat(arrayfun(@(x)[x;x],timesSol.steps,'UniformOutput',false));
        timesSol.tPieceSteps = timesSol.tPieceSteps(2:end-1);
        
    
    end
    
    for ci = 1:numel(wellSols{1})
        if opt.wc
            figure(figN); figN = figN+1;
            subplot(3,1,1)
            qls = qOs(ci,:)+qWs(ci,:);
            if simulateFlag
                qlsS = qOsS(ci,:)+qWsS(ci,:);
                plot(times.tPieceSteps, qls, 'bx-',timesSol.tPieceSteps, qlsS, 'ro-')
                legend('MS','Fwd')
            else
                plot(times.tPieceSteps, qls, 'x-')
            end
            ylabel('q_l (m^3/day)');
            title(strcat('Well: ',wellSols{1}(ci).name,' Type: ',wellSols{1}(ci).type,' Sign: ',intType2stringType(wellSols{1}(ci).sign)) );
            
            subplot(3,1,2)
            wcuts = qWs(ci,:)./qls;
            if simulateFlag
                wcutsS = qWsS(ci,:)./qlsS;
                plot(times.tPieceSteps, wcuts, 'bx-',timesSol.tPieceSteps, wcutsS, 'ro-')
                legend('MS','Fwd')
            else
                plot(times.tPieceSteps, wcuts, 'x-')
            end
            ylabel('WCUT');
            
            
        else
            figure(figN); figN = figN+1;
            subplot(3,1,1)
            if simulateFlag
                plot(times.tPieceSteps, qOs(ci,:), 'bx-',timesSol.tPieceSteps, qOsS(ci,:), 'ro-')
                legend('MS','Fwd')
            else
                plot(times.tPieceSteps, qOs(ci,:), 'x-')
            end
            ylabel('q_o (m^3/day)');
            title(strcat('Well: ',wellSols{1}(ci).name,' Control target: ',wellSols{1}(ci).type,' Type: ',intType2stringType(wellSols{1}(ci).sign)) );
            
            subplot(3,1,2)
            if simulateFlag
                plot(times.tPieceSteps, qWs(ci,:), 'bx-',timesSol.tPieceSteps, qWsS(ci,:), 'ro-')
                legend('MS','Fwd')
            else
                plot(times.tPieceSteps, qWs(ci,:), 'x-')
            end
            ylabel('q_w (m^3/day)');
            
            
        end
        subplot(3,1,3)
        if simulateFlag
            plot(times.tPieceSteps, bhp(ci,:), 'bx-',timesSol.tPieceSteps, bhpS(ci,:), 'ro-')
            legend('MS','Fwd')
        else
            plot(times.tPieceSteps, bhp(ci,:), 'x-')
        end
        ylabel('bhp (bar)');
        xlabel('time(day)')
        
        
    end
end
    
% plotting network constraints
if opt.plotNetsol
    numNetworkConstraints = 5;
    figN = plotNetworkConstraints(v, lbv, ubv, nScale, times, numNetworkConstraints, figN);    
end

if opt.plotNetControls
    numControls = 1;   
    
    figN = plotNetworkControls(u, lbu, ubu, uScalePlot, times, numControls,  figN);
end

if opt.plotObjective
    %if opt.plotObj
totalPredictionSteps = getTotalPredictionSteps(ss);
    objs = zeros(1,totalPredictionSteps);
    for k = 1:totalPredictionSteps
        
        if k > 1
            [objs(k)] = callArroba(obj{k},{x{k-1},u{callArroba(ss.ci,{k})},v{k}});
        elseif k == 1
            [objs(k)] = callArroba(obj{k},{ss.state,u{callArroba(ss.ci,{k})},v{k}});
        else
            error('what?')
        end
    end
    totalObj = -sum(objs);
    objDay =  arrayfun(@(x,y)-x/y,objs,diff(times.steps)');
    
    
    %plot results
    
    objDay = cell2mat(arrayfun(@(x)[x,x],objDay,'UniformOutput',false));
    
    figure(figN); figN = figN+1;
    plot(times.tPieceSteps, objDay, '-x')
    xlabel('time (day)');
    title(strcat('Objective/day. Total Objective: ',num2str(totalObj)) );
end




end
