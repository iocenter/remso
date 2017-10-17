function [  ] = plotSolutionOW( x,u,v,d, lbv, ubv, lbu, ubu, ss,obj,times,xScale,uScale,vScale, nScale, uScalePlot,schedules,wellSol, netSol, lbuPot,ubuPlot,ulim,minState,maxState,varargin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% varargin = {'simulate',[],'xScale',xScale,'uScale',cellControlScales,'uScalePlot',cellControlScalesPlot,'schedules',mShootingP.schedules}
% opt = struct('simulate',[],'simFlag',false,'plotWellSols',true,'plotSchedules',true,'plotObjective',true,'pF',@(x)x,'sF',@(x)x,'figN',1000,'wc',false,'reservoirP',[],'plotSweep',false);
opt = struct('simulate',[],'simFlag',false,'plotWellSols',true, 'plotNetsol', false, 'numNetConstraints', 0, 'plotNetControls', false, 'numNetControls', 0, 'plotCumulativeObjective', false, ...
            'plotSchedules',true,'plotObjective',true,'pF',@(x)x,'sF',@(x)x,'figN',1000,'wc',false,'reservoirP',[],'plotSweep',false, ...
            'freqCst', 0, 'pressureCst', 0, 'flowCst', 0, 'fixedWells', [], 'stepBreak', numel(v), 'extremePoints', [], ...
            'qlMin', [], 'qlMax', [], 'nStages', [], 'freqMin', [], 'freqMax', [], 'baseFreq', [], ...
            'plotNetwork', false, 'dpFunction', [], 'plotChokeConstraints', false);
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

nfw = numel(opt.fixedWells);
nwc  = cellfun(@(s) s - nfw, nW,'UniformOutput',false);

w = cellfun(@(ui,nw)ui(1:nw),u,nwc,'UniformOutput',false);
p = cellfun(@(ui,nw)ui(nw+1:end),u,nwc,'UniformOutput',false);

wScale = cellfun(@(ui,nw)ui(1:nw),uScale,nwc,'UniformOutput',false);
wScalePlot = cellfun(@(ui,nw)ui(1:nw),uScalePlot,nW,'UniformOutput',false);

pScale = cellfun(@(ui,nw)ui(nw+1:end),uScale,nW,'UniformOutput',false);

if opt.plotSchedules
    [uM,schedulesSI] = scaleSchedulePlot(w,schedules,wScale,wScalePlot, 'fixedWells' ,opt.fixedWells);
    
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
    [uM,schedulesSI] = scaleSchedulePlot(w,schedules,wScale,wScalePlot, 'fixedWells', opt.fixedWells);
    
    
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
    
    
    if opt.plotNetwork
        netSols= cell(totalPredictionSteps,1);
        shootingGuess = cell(totalPredictionSteps,1);
        if ~isempty(v)
            for i=1:totalPredictionSteps
                [shootingGuess{i}.wellSol] = algVar2wellSol( v{i},wellSol,'vScale', vScale,...
                    'activeComponents',opt.reservoirP.system.activeComponents);

                netSols{i} = runNetwork(netSol, shootingGuess{i}.wellSol, [],  'fluid',opt.reservoirP.fluid, 'activeComponents',opt.reservoirP.system.activeComponents, 'dpFunction', opt.dpFunction);
            end
        end          
          
        qlMin = cellfun(@(x) opt.qlMin, netSols, 'UniformOutput', false);
        qlMax = cellfun(@(x) opt.qlMax, netSols, 'UniformOutput', false);
        numStages = cellfun(@(x) opt.nStages, netSols, 'UniformOutput', false);
        fref = cellfun(@(x) opt.baseFreq, netSols, 'UniformOutput', false);
        
        [freq, qf,  qpump_min, qpump_max,dhf, ~] = cellfun(@nonlinearPumpConstraints, netSols,  fref, numStages, qlMin, qlMax,  'UniformOutput', false);
        pressures = cell2mat((cellfun(@(ns) ns.pV, netSols, 'UniformOutput', false))');
        eqp = getEdge(netSol, netSol.Eeqp);
        time = cumsum(opt.reservoirP.schedule.step.val)./day; 
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
            waterFlow = qWs(ci,:);
            if simulateFlag
                plot(times.tPieceSteps, waterFlow, 'bx-',timesSol.tPieceSteps, qWsS(ci,:), 'ro-')
                legend('MS','Fwd')
            else
                plot(times.tPieceSteps, waterFlow, 'x-')
            end
            ylabel('q_w (m^3/day)');
            if (max(waterFlow) - min(waterFlow)) < 1e-2 *meter^3/day %% all elements are equal
                prevAx = axis;                
                axis([prevAx(1) prevAx(2) max(waterFlow(1)-1,0) waterFlow(1)+1]);
            end
            
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
        
        if opt.plotNetwork
            if any(netSol.Vsrc==ci)              
                figure(figN); figN = figN+1;     

                vi = getVertex(netSol, ci);
                branch = [vi.id];
                while ~isempty(vi.Eout)
                    ei = getEdge(netSol, vi.Eout);
                    vi = getVertex(netSol, ei.vout);
                    branch = [branch; vi.id];
                end

                for k=1:numel(branch)-1                    
                    if any(branch(k) == vertcat(eqp.vin))
                        pumpInd = [branch(k); branch(k+1)];
                        if pressures(branch(k),:) > pressures(branch(k+1),:)
                            warning('pump decreasing pressure');
                        end
                    else
                        if pressures(branch(k),:) < pressures(branch(k+1),:)
                            warning('pipeline increasing pressure');
                        end
                    end
                end                
                
                subplot(2,1, 1);
                freqPump =  cellfun(@(x) x(ci-numel(netSol.VwInj)), freq); 
                plot(time, opt.freqMin(ci-numel(netSol.VwInj)), 'rv-', time, freqPump, '-x', time, opt.freqMax(ci-numel(netSol.VwInj)), 'r^-');
                xlabel('time (day)');
                ylabel('frequency (Hz)');
                title(strcat('Pump Frequency: ', ' ',  netSol.V(ci).name));
                ylim([0 max(opt.freqMax(ci-numel(netSol.VwInj)))+10])                
                
                
                subplot(2,1,2);                
                qliqPump =  cellfun(@(x) x(ci-numel(netSol.VwInj)), qf);
                qlMin    =  cellfun(@(x) x(ci-numel(netSol.VwInj)), qpump_min);
                qlMax    =  cellfun(@(x) x(ci-numel(netSol.VwInj)), qpump_max);                
                
                plot(time,qlMin./(meter^3/day) ,'rv-', time, -qliqPump./(meter^3/day), '-x', time, qlMax./(meter^3/day) ,'rv-');
                xlabel('time (day)');
                ylabel('liq. flow  rate (sm3/d)');
                title(strcat('Pump Liq. Flow Rate: ', ' ',  netSol.V(ci).name));
                
                
                figure(figN); figN = figN+1;                     
                plot(time, pressures(branch, :)./barsa, '-x');
                xlabel('time (day)');
                ylabel('pressure (bar)');
                title(strcat('Network Branch of Well: ', ' ',  netSol.V(ci).name));
                
                legend('p(v1) - Source','p (v2) - Ups. Equip.', ' p (v3) - Downs. Equip.', 'p (v4) - Well-Head', 'p (v5) - Manifold', 'p (v6) -- Sink');                
                
                figure(figN); figN = figN+1; 
                wellId =ci-numel(netSol.VwInj);
                numFreq = 60;
                numFlows = 55;                
                f0 = 60;                
                freqStep = (opt.freqMax(wellId)-opt.freqMin(wellId))/numFreq;
                qmin = opt.qlMin(wellId)./(meter^3/day);
                qmax = opt.qlMax(wellId)./(meter^3/day);
                
                
                dhfk     = zeros(numFreq,numFlows);                
                
                i=0;                   
                for fi=opt.freqMin(wellId):freqStep:opt.freqMax
                    i = i+1;
                    k=0; 
                    qminF = qmin*(fi/f0);
                    qmaxF = qmax*(fi/f0);   
                    flowStep = (qmaxF-qminF)/numFlows;
                    for qk=qminF:flowStep:qmaxF
                        k = k+1;
                        
                        a0 = 19.37;       % m
                        a1 = -1.23e-02;     % in m/sm3/d
                        a2 = 2.24e-05;      % in m/(sm3/d)^2
                        a3 = -1.86e-08;     % in m/(sm3/d)^3
                        a4 = 4.13e-12;      % in m/(sm3/d)^4
                        
                        dhf0 = a4*(qk*(f0/fi)).^4 + a3*(qk*(f0/fi))^3 + a2*(qk*(f0/fi))^2 + a1*(qk*(f0/fi)) + a0;
                        dhfk(i,k) = dhf0*((fi/f0)^2)*opt.nStages(wellId);   
                    end                    
                end    
                
                i=0;
                qminF = zeros(numel(numFreq),1);
                qmaxF = zeros(numel(numFreq),1);
                for fi=opt.freqMin(wellId):freqStep:opt.freqMax
                    i = i+1;
                    qminF(i) = qmin*(fi/f0);
                    qmaxF(i) = qmax*(fi/f0);             
                end
                line(qminF, dhfk(:,1), 'LineStyle','-', 'Color',[0 0 0], 'LineW', 1.5);
                hold on;
                
                line(qmaxF, dhfk(:,end), 'LineStyle','-', 'Color',[0 0 0], 'LineW', 1.5);
                hold on;
               
                for fi=[opt.freqMin(wellId), opt.freqMax(wellId)]                     
                     qminF = qmin*(fi/f0);
                     qmaxF = qmax*(fi/f0);   
                     flowStep = (qmaxF-qminF)/numFlows;
                     if fi==opt.freqMin(wellId)
                        plot(qminF:flowStep:qmaxF, dhfk(1,:), 'Color',[0 0 0], 'LineW', 1.5);
                     else
                         plot(qminF:flowStep:qmaxF, dhfk(end,:), 'Color',[0 0 0], 'LineW', 1.5);
                     end
                end
                
                liqPump = (-qliqPump./(meter^3/day))';
                dhfPump = cellfun(@(x) x(ci-numel(netSol.VwInj)), dhf);
                dhPump = dhfPump';
                zPump = zeros(size(liqPump));
                
                col = time';  % This is the color, vary with x in this case.                
                
                surface([liqPump;liqPump],[dhPump;dhPump],[zPump;zPump],[col;col],...
                    'facecol','no',...
                    'edgecol','interp',...
                    'linew',1.5);
                c = colorbar;
                c.Label.String = 'time (days)';
                
                hold all;
                
                xlabel('local liq. flow rate (sm3/d)');
                ylabel('pump head (m)');
                title(strcat(netSol.V(ci).name, ': ESP Operating Envelope.'));
            end
        end
       
    end
end
    
% plotting network constraints
if opt.plotNetsol
    numNetworkConstraints = opt.numNetConstraints;
    if opt.plotChokeConstraints
        n = cellfun(@(vk) vk(end-numNetworkConstraints:end-1), v, 'UniformOutput', false);
        lbn  = cellfun(@(vk) vk(end-numNetworkConstraints:end-1), lbv, 'UniformOutput', false);     
        ubn  = cellfun(@(vk) vk(end-numNetworkConstraints:end-1), ubv, 'UniformOutput', false);
        
        dp = zeros(numNetworkConstraints, numel(v));
        ldp = zeros(numNetworkConstraints, numel(v));
        udp = zeros(numNetworkConstraints, numel(v));
        
        for i=1:numel(v)
            dp(:,i)  = n{i}.*nScale;
            ldp(:,i) = lbn{i}.*nScale;
            udp(:,i) = ubn{i}.*nScale;            
        end
        
       figN = plotPressureConstraints(dp, ldp, udp, times, numNetworkConstraints, figN);
    else
        figN = plotNetworkConstraints(v, lbv, ubv, nScale, times, numNetworkConstraints, figN, ...
            'freqCst', opt.freqCst, ...
            'pressureCst', opt.pressureCst, ...
            'flowCst', opt.flowCst, ...
            'nW', cell2mat(nW), ...
            'vScale', vScale, ...
            'extremePoints', opt.extremePoints, ...
            'qlMin', opt.qlMin, ...
            'qlMax', opt.qlMax, ...
            'nStages', opt.nStages, ...
            'freqMin', opt.freqMin, ...
            'freqMax', opt.freqMax, ...
            'baseFreq',opt.baseFreq);
    end
end

if opt.plotNetControls
    numControls = opt.numNetControls;   
    
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

if opt.plotCumulativeObjective    
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
    
    %plot results
    
    figure(figN); 
    figN = figN+1;
    plot(times.steps(2:end), cumsum(-objs), '-x');
%     plot(times.tPieceSteps, objAcum, '-x')
    xlabel('time (day)');
    title(strcat('Cumulative Objective. Total Objective: ',num2str(totalObj)) );
end

end
