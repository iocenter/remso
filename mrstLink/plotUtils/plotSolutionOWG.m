function [  ] = plotSolutionOWG( x,u,v,d,ss,obj,times,xScale,uScale,vScale,uScalePlot,schedules,wellSol,lbuPot,ubuPlot,ulim,minState,maxState,varargin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% varargin = {'simulate',[],'xScale',xScale,'uScale',cellControlScales,'uScalePlot',cellControlScalesPlot,'schedules',mShootingP.schedules}
opt = struct('simulate',[],'simFlag',false,'plotWellSols',true,'plotSchedules',true,'plotObjective',true,'pF',@(x)x,'sF',@(x)x,'figN',1000,'wc',false,...
    'activeComponents',struct('oil',1,'water',1,'gas',1,'polymer',0,'disgas',0,'vapoil',0,'T',0,'MI',0),...% default OWG
    'fluid',[],...
    'system',[]);
opt = merge_options(opt, varargin{:});

figN = opt.figN;

disgas = opt.activeComponents.disgas;
vapoil = opt.activeComponents.vapoil;



xM = cellfun(@(xi)stateVector2stateMrst( xi,'xScale',xScale,...
    'activeComponents',opt.activeComponents,...
    'fluid',opt.fluid,...
    'system',opt.system),...
    [ss.state;x],'UniformOutput',false);


nx = numel(x{1})/3;

[p,sW,rGH] = cellfun(@(xi)stateVector2PsWrGH(xi,xScale,nx),d,'UniformOutput',false);


pPlot  = cellfun(@(x)opt.pF(x.pressure/barsa),xM,'UniformOutput',false)';
dpPlot = cellfun(@(x)opt.pF(x/barsa),p,'UniformOutput',false)';

sWPlot  = cellfun(@(x)opt.sF(x.s(:,1)),xM,'UniformOutput',false)';
dsWPlot = cellfun(@(x)opt.sF(x),sW,'UniformOutput',false)';

sGPlot  = cellfun(@(x)opt.sF(x.s(:,3)),xM,'UniformOutput',false)';

if disgas
    rsPlot  = cellfun(@(x)opt.sF(x.rs),xM,'UniformOutput',false)';
end
if vapoil
    rvPlot  = cellfun(@(x)opt.sF(x.rv),xM,'UniformOutput',false)';
end


drGHPlot = cellfun(@(x)opt.sF(x),rGH,'UniformOutput',false)';



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
plot(times.steps,cell2mat(sWPlot),'-x');
ylabel('Water saturation')
xlabel('time (days)')
ylim([minState.sW,maxState.sW])
subplot(2,1,2)
plot(times.steps(2:end),cell2mat(dsWPlot),'-x');
ylabel('Water saturation error')
xlabel('time (days)')



figure(figN); figN = figN+1;
subplot(2,1,1)
plot(times.steps,cell2mat(sGPlot),'-x');
ylabel('Gas saturation')
xlabel('time (days)')
ylim([0,1]);
subplot(2,1,2)
plot(times.steps(2:end),cell2mat(drGHPlot),'-x');
ylabel('rGH  error')
xlabel('time (days)')

if disgas || vapoil
    figure(figN); figN = figN+1;
end
if disgas
    subplot(disgas+vapoil,1,1)
    plot(times.steps,cell2mat(rsPlot),'-x');
    ylabel('rs')
    xlabel('time (days)')
    rsMax = opt.fluid.rsSat(maxState.pressure);
    ylim([0,rsMax])
end
if vapoil
    subplot(disgas+vapoil,1,disgas+vapoil)
    plot(times.steps,cell2mat(rvPlot),'-x');
    ylabel('rv')
    xlabel('time (days)')
    rvMax = opt.fluid.rvSat(maxState.pressure);
    ylim([0,rvMax])
end



if opt.plotSchedules
    [uM,schedulesSI] = scaleSchedulePlot(u,schedules,uScale,uScalePlot);
    
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
    [uM,schedulesSI] = scaleSchedulePlot(u,schedules,uScale,uScalePlot);
    
    
    wellSols = cellfun(@(vk)algVar2wellSol(vk,wellSol,'vScale',vScale,...
                                           'activeComponents',opt.activeComponents),...
                       v,'UniformOutput',false )';
    
    [qWs, qOs, qGs, bhp] = wellSolToVector(wellSols);
    qWs = cell2mat(arrayfun(@(x)[x,x],qWs'*day,'UniformOutput',false));
    qOs = cell2mat(arrayfun(@(x)[x,x],qOs'*day,'UniformOutput',false));
    qGs = cell2mat(arrayfun(@(x)[x,x],qGs'*day,'UniformOutput',false));
    bhp = cell2mat(arrayfun(@(x)[x,x],bhp'/barsa,'UniformOutput',false));
    
    simulateFlag = ~isempty(opt.simulate) && opt.simFlag;
    if simulateFlag
        schedule = mergeSchedules(schedulesSI);
        wellSolsS = opt.simulate(schedule);
        [qWsS, qOsS, qGsS, bhpS] = wellSolToVector(wellSolsS);
        qWsS = cell2mat(arrayfun(@(x)[x,x],qWsS'*day,'UniformOutput',false));
        qOsS = cell2mat(arrayfun(@(x)[x,x],qOsS'*day,'UniformOutput',false));
        qGsS = cell2mat(arrayfun(@(x)[x,x],qGsS'*day,'UniformOutput',false));
        bhpS = cell2mat(arrayfun(@(x)[x,x],bhpS'/barsa,'UniformOutput',false));
    end
    
    for ci = 1:size(wellSols{1},2)
        figure(figN); figN = figN+1;
        if opt.wc
            subplot(4,1,1)
            qls = qOs(ci,:)+qWs(ci,:);
            if simulateFlag
                qlsS = qOsS(ci,:)+qWsS(ci,:);
                plot(times.tPieceSteps, qls, 'bx-',times.tPieceSteps, qlsS, 'ro-')
                legend('MS','Fwd')
            else
                plot(times.tPieceSteps, qls, 'x-')
            end
            ylabel('q_l (meter^3/day)');
            title(strcat('Well: ',wellSols{1}(ci).name,' Type: ',wellSols{1}(ci).type,' Sign: ',intType2stringType(wellSols{1}(ci).sign)) );
            
            subplot(4,1,2)
            wcuts = qWs(ci,:)./qls;
            if simulateFlag
                wcutsS = qWsS(ci,:)./qlsS;
                plot(times.tPieceSteps, wcuts, 'bx-',times.tPieceSteps, wcutsS, 'ro-')
                legend('MS','Fwd')
            else
                plot(times.tPieceSteps, wcuts, 'x-')
            end
            ylabel('WCUT');
            
            
        else
            subplot(4,1,1)
            if simulateFlag
                plot(times.tPieceSteps, qOs(ci,:), 'bx-',times.tPieceSteps, qOsS(ci,:), 'ro-')
                legend('MS','Fwd')
            else
                plot(times.tPieceSteps, qOs(ci,:), 'x-')
            end
            ylabel('q_o (meter^3/day)');
            title(strcat('Well: ',wellSols{1}(ci).name,' Type: ',wellSols{1}(ci).type,' Sign: ',intType2stringType(wellSols{1}(ci).sign)) );
            
            subplot(4,1,2)
            if simulateFlag
                plot(times.tPieceSteps, qWs(ci,:), 'bx-',times.tPieceSteps, qWsS(ci,:), 'ro-')
                legend('MS','Fwd')
            else
                plot(times.tPieceSteps, qWs(ci,:), 'x-')
            end
            ylabel('q_w (meter^3/day)');
            
            
        end
        
        subplot(4,1,3)
        if simulateFlag
            plot(times.tPieceSteps, qGs(ci,:), 'bx-',times.tPieceSteps, qGsS(ci,:), 'ro-')
            legend('MS','Fwd')
        else
            plot(times.tPieceSteps, qGs(ci,:), 'x-')
        end
        ylabel('q_g (meter^3/day)');
        xlabel('time(day)')
        
        
        subplot(4,1,4)
        if simulateFlag
            plot(times.tPieceSteps, bhp(ci,:), 'bx-',times.tPieceSteps, bhpS(ci,:), 'ro-')
            legend('MS','Fwd')
        else
            plot(times.tPieceSteps, bhp(ci,:), 'x-')
        end
        ylabel('bhp (barsa)');
        xlabel('time(day)')
        
    end
else
    
    
    
    
end

if opt.plotObjective
    %if opt.plotObj
    totalPredictionSteps = getTotalPredictionSteps(ss);
    objs = zeros(1,totalPredictionSteps);
    for k = 1:totalPredictionSteps
        
        if k > 1
            [objs(k)] = callArroba(obj{k},{x{k-1},u{ss.ci(k)},v{k}});
        elseif k == 1
            [objs(k)] = callArroba(obj{k},{ss.state,u{ss.ci(k)},v{k}});
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



function [p,sW,rGH] = stateVector2PsWrGH(stateVector,xScale,nx)

if ~isempty(xScale)
    stateVector = stateVector.*xScale;
end



p = stateVector(1:nx);
sW = stateVector(nx+1:2*nx);
rGH = stateVector(2*nx+1:end);


end