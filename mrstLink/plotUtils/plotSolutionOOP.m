function [  ] = plotSolutionOOP( x,u,v,xd,ss,obj,model,uScale,stepSchedules,controlSchedules,minState,maxState,varargin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

opt = struct('simFlag',false,'plotWellSols',true,'plotStateErros',true,'plotObjective',true,'pF',@(x)x,'sF',@(x)x,'figN',1000,'wc',false,'units','METRIC');
opt = merge_options(opt, varargin{:});

totalPredictionSteps = numel(x);

figN = opt.figN;

[units,unitsName] = unit_system(opt.units);

times.steps = [stepSchedules(1).time;arrayfun(@(x)(x.time+sum(x.step.val)),stepSchedules)]/units.time;
times.tPieceSteps = cell2mat(arrayfun(@(x)[x;x],times.steps,'UniformOutput',false));
times.tPieceSteps = times.tPieceSteps(2:end-1);

times.controls = [controlSchedules(1).time;arrayfun(@(x)(x.time+sum(x.step.val)),controlSchedules)]/units.time;
times.tPieceControls = cell2mat(arrayfun(@(x)[x;x],times.controls,'UniformOutput',false));
times.tPieceControls = times.tPieceControls(2:end-1);


if model.gas
    disgas = model.disgas;
    vapoil = model.vapoil;
else
    disgas = false;
    vapoil = false;
end

xM = cellfun(@(xi)model.toMRSTStates(xi),[ss.state;x],'UniformOutput',false);


pPlot  = cellfun(@(x)opt.pF(model.getProp(x,'pressure')/units.press),xM,'UniformOutput',false)';
sWPlot  = cellfun(@(x)opt.sF(model.getProp(x,'sw')),xM,'UniformOutput',false)';
sGPlot  = cellfun(@(x)opt.sF(model.getProp(x,'sg')),xM,'UniformOutput',false)';

if disgas
    rsPlot  = cellfun(@(x)opt.sF(model.getProp(x,'rs')),xM,'UniformOutput',false)';
end
if vapoil
    rvPlot  = cellfun(@(x)opt.sF(model.getProp(x,'rv')),xM,'UniformOutput',false)';
end


figure(figN); figN = figN+1;
plot(times.steps,cell2mat(pPlot),'-x');
ylabel(['Pressure (' unitsName.press ')'])
xlabel(['time (',unitsName.time,')'])
ylim([minState.pressure,maxState.pressure]/units.press)


figure(figN); figN = figN+1;
plot(times.steps,cell2mat(sWPlot),'-x');
ylabel('Water saturation')
xlabel(['time (',unitsName.time,')'])
ylim([minState.sW,maxState.sW])
xlabel(['time (',unitsName.time,')'])


if model.gas
    figure(figN); figN = figN+1;
    plot(times.steps,cell2mat(sGPlot),'-x');
    ylabel('Gas saturation')
    xlabel(['time (',unitsName.time,')'])
    ylim([0,1]);
end


if disgas
    figure(figN); figN = figN+1;

    plot(times.steps,cell2mat(rsPlot),'-x');
    ylabel('rs')
    xlabel(['time (',unitsName.time,')'])
    rsMax = opt.fluid.rsSat(maxState.pressure);
    ylim([0,rsMax])
end
if vapoil
	figure(figN); figN = figN+1;

    plot(times.steps,cell2mat(rvPlot),'-x');
    ylabel('rv')
    xlabel(['time (',unitsName.time,')'])
    rvMax = model.fluid.rvSat(maxState.pressure);
    ylim([0,rvMax])
end

if opt.plotStateErros
    nx = model.G.cells.num;

    if model.gas 
        [x1,x2,x3] = cellfun(@(xi)splitStateVector(xi,nx,model.gas),xd,'UniformOutput',false);
    else
        [x1,x2] = cellfun(@(xi)splitStateVector(xi,nx,model.gas),xd,'UniformOutput',false);
    end
    
 	figure(figN); figN = figN+1;
    plot(times.steps(2:end),cell2mat(x1'),'-x');
    ylabel('Scaled state 1')
    xlabel(['time (',unitsName.time,')'])
    
 	figure(figN); figN = figN+1;
    plot(times.steps(2:end),cell2mat(x2'),'-x');
    ylabel('Scaled state 2')
    xlabel(['time (',unitsName.time,')'])
    
    if model.gas
        figure(figN); figN = figN+1;
        plot(times.steps(2:end),cell2mat(x3'),'-x');
        ylabel('Scaled state 3')
        xlabel(['time (',unitsName.time,')']) 
    end

    
    
end


if opt.plotWellSols
    
    usliced = cell(totalPredictionSteps,1);
    uScaleSliced = cell(totalPredictionSteps,1);
    for k = 1:totalPredictionSteps
        usliced{k} = u{callArroba(ss.ci,{k})};
        uScaleSliced{k} = uScale{callArroba(ss.ci,{k})};
    end    
    
    schedules = mat2cell(stepSchedules,ones(numel(stepSchedules),1),1);
    schedules = cellfun(@(ui,si,usi)u2schedules( {ui},si,'uScale',{usi}),usliced,schedules,uScaleSliced,'UniformOutput',false);
    
    lastWells = cellfun(@(s)s.control(end).W,schedules,'UniformOutput',false);
    
    
    wellSols = cellfun(@(vk,w)model.toWellSol(vk,w),...
                       v,lastWells,'UniformOutput',false )';
                   
    if ~model.gas
        wellSols = cellfun(@(x)arrayfun(@(y)subsasgn(y,struct('type','.','subs','qGs'),0),x),wellSols,'UniformOutput',false);
    end                   
                   
    
    [qWs, qOs, qGs, bhp] = wellSolToVector(wellSols);
    qWs = cell2mat(arrayfun(@(x)[x,x],qWs'/(units.liqvol_s/units.time),'UniformOutput',false));
    qOs = cell2mat(arrayfun(@(x)[x,x],qOs'/(units.liqvol_s/units.time),'UniformOutput',false));
    qGs = cell2mat(arrayfun(@(x)[x,x],qGs'/(units.gasvol_s/units.time),'UniformOutput',false));
    bhp = cell2mat(arrayfun(@(x)[x,x],bhp'/(units.press),'UniformOutput',false));
    
    simulateFlag = opt.simFlag;
    if simulateFlag
        schedule = mergeSchedules([schedules{:}]);            
        initState = model.toMRSTStates(ss.state);
        
        [wellSolsS,~,schedulereport] = simulateScheduleAD(initState, model, schedule, 'OutputMinisteps',true);
        [scheduleS] = convertReportToSchedule(schedulereport,schedule);

        if ~model.gas
            wellSolsS = cellfun(@(x)arrayfun(@(y)subsasgn(y,struct('type','.','subs','qGs'),0),x),wellSolsS,'UniformOutput',false);
        end       

        [qWsS, qOsS, qGsS, bhpS] = wellSolToVector(wellSolsS);
        qWsS = cell2mat(arrayfun(@(x)[x,x],qWsS'/(units.liqvol_s/units.time),'UniformOutput',false));
        qOsS = cell2mat(arrayfun(@(x)[x,x],qOsS'/(units.liqvol_s/units.time),'UniformOutput',false));
        qGsS = cell2mat(arrayfun(@(x)[x,x],qGsS'/(units.gasvol_s/units.time),'UniformOutput',false));
        bhpS = cell2mat(arrayfun(@(x)[x,x],bhpS'/(units.press),'UniformOutput',false));
        
        steps = [scheduleS.time,scheduleS.time+cumsum(scheduleS.step.val)']/(units.time)  ;      
        tPieceSteps = cell2mat(arrayfun(@(x)[x;x],steps,'UniformOutput',false));
        tPieceSteps = tPieceSteps(2:end-1);        
    end
    
    g = model.gas;
    nPlots = 3+g;
    for ci = 1:size(wellSols{1},2)
        figure(figN); figN = figN+1;
        if opt.wc
            subplot(nPlots,1,1)
            qls = qOs(ci,:)+qWs(ci,:);
            if simulateFlag
                qlsS = qOsS(ci,:)+qWsS(ci,:);
                plot(times.tPieceSteps, qls, 'bx-',tPieceSteps, qlsS, 'ro-')
                legend('MS','Fwd')
            else
                plot(times.tPieceSteps, qls, 'x-')
            end
            ylabel(['q_l (' unitsName.liqvol_s '/' unitsName.time  ')']);
            title(strcat('Well: ',wellSols{1}(ci).name,' Type: ',wellSols{1}(ci).type,' Sign: ',intType2stringType(wellSols{1}(ci).sign)) );
            
            subplot(nPlots,1,2)
            wcuts = qWs(ci,:)./qls;
            if simulateFlag
                wcutsS = qWsS(ci,:)./qlsS;
                plot(times.tPieceSteps, wcuts, 'bx-',tPieceSteps, wcutsS, 'ro-')
                legend('MS','Fwd')
            else
                plot(times.tPieceSteps, wcuts, 'x-')
            end
            ylabel('WCUT');
            
            
        else
            subplot(nPlots,1,1)
            if simulateFlag
                plot(times.tPieceSteps, qOs(ci,:), 'bx-',tPieceSteps, qOsS(ci,:), 'ro-')
                legend('MS','Fwd')
            else
                plot(times.tPieceSteps, qOs(ci,:), 'x-')
            end
            ylabel(['q_o (' unitsName.liqvol_s '/' unitsName.time  ')']);
            title(strcat('Well: ',wellSols{1}(ci).name,' Type: ',wellSols{1}(ci).type,' Sign: ',intType2stringType(wellSols{1}(ci).sign)) );
            
            subplot(nPlots,1,2)
            if simulateFlag
                plot(times.tPieceSteps, qWs(ci,:), 'bx-',tPieceSteps, qWsS(ci,:), 'ro-')
                legend('MS','Fwd')
            else
                plot(times.tPieceSteps, qWs(ci,:), 'x-')
            end
            ylabel(['q_w (' unitsName.liqvol_s '/' unitsName.time  ')']);
            
            
        end
        
        if model.gas
            subplot(nPlots,1,3)
            if simulateFlag
                plot(times.tPieceSteps, qGs(ci,:), 'bx-',tPieceSteps, qGsS(ci,:), 'ro-')
                legend('MS','Fwd')
            else
                plot(times.tPieceSteps, qGs(ci,:), 'x-')
            end
            ylabel(['q_g (' unitsName.liqvol_s '/' unitsName.time  ')']);
            xlabel(['time('  unitsName.time  ')'])
        end
        
        
        subplot(nPlots,1,nPlots)
        if simulateFlag
            plot(times.tPieceSteps, bhp(ci,:), 'bx-',tPieceSteps, bhpS(ci,:), 'ro-')
            legend('MS','Fwd')
        else
            plot(times.tPieceSteps, bhp(ci,:), 'x-')
        end
        ylabel(['bhp (' unitsName.press ')']);
        xlabel(['time('  unitsName.time  ')'])
        
    end
else
    
    
    
    
end

if opt.plotObjective
    %if opt.plotObj
    totalPredictionSteps = getTotalPredictionSteps(ss);
    objs = zeros(1,totalPredictionSteps);
    for k = 1:totalPredictionSteps
        cik = callArroba(ss.ci,{k});
        if k > 1
            [objs(k)] = callArroba(obj{k},{x{k-1},u{cik},v{k}});
        elseif k == 1
            [objs(k)] = callArroba(obj{k},{ss.state,u{cik},v{k}});
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
    xlabel(['time ('  unitsName.time  ')']);
    title(strcat(['Objective/' unitsName.time '. Total Objective: '],num2str(totalObj)) );
end




end



function [p,sW,rGH] = splitStateVector(stateVector,nx,gas)

p = stateVector(1:nx);
sW = stateVector(nx+1:2*nx);
if gas    
    rGH = stateVector(2*nx+1:end);
else
    rGH = [];
end


end

%% copied from mrst/modules/deckformat/deckinput/convertDeckUnits.m
function [u,n] = unit_system(rspec)
   metric = strcmp(rspec, 'METRIC');
   field  = strcmp(rspec, 'FIELD');
   lab    = strcmp(rspec, 'LAB');
   SI     = strcmp(rspec, 'SI');

   if sum([metric, field, lab, SI]) ~= 1,
      error(id('USys:Unknown'), ...
            'Input unit system must be either METRIC, FIELD, LAB, or SI.');
   end

   if metric,
      u = struct('length'   , meter             , ...
                 'time'     , day               , ...
                 'density'  , kilogram / meter^3, ...
                 'press'    , barsa             , ...
                 'concentr' , kilogram / meter^3, ... % Concentration
                 'compr'    , 1 / barsa         , ... % Compressibility
                 'viscosity', centi*poise       , ...
                 'perm'     , milli*darcy       , ...
                 'liqvol_s' , meter^3           , ... % Liquid vol, surf
                 'liqvol_r' , meter^3           , ... % Liquid vol, res
                 'gasvol_s' , meter^3           , ... % Gas vol, surf
                 'gasvol_r' , meter^3           , ... % Gas vol, res
                 'trans'    , centi*poise * meter^3 / (day * barsa));
      n = struct('length'   , 'meter'             , ...
                 'time'     , 'day'               , ...
                 'density'  , 'kilogram / meter^3', ...
                 'press'    , 'barsa'             , ...
                 'concentr' , 'kilogram / meter^3', ... % Concentration
                 'compr'    , '1 / barsa'         , ... % Compressibility
                 'viscosity', 'centi*poise'       , ...
                 'perm'     , 'milli*darcy'       , ...
                 'liqvol_s' , 'meter^3'           , ... % Liquid vol, surf
                 'liqvol_r' , 'meter^3'           , ... % Liquid vol, res
                 'gasvol_s' , 'meter^3'           , ... % Gas vol, surf
                 'gasvol_r' , 'meter^3'           , ... % Gas vol, res
                 'trans'    , 'centi*poise * meter^3 / (day * barsa)');
   elseif field
      u = struct('length'   , ft                , ...
                 'time'     , day               , ...
                 'density'  , pound / ft^3      , ...
                 'press'    , psia              , ...
                 'concentr' , pound / stb       , ... % Concentration
                 'compr'    , 1 / psia          , ...
                 'viscosity', centi*poise       , ...
                 'perm'     , milli*darcy       , ...
                 'liqvol_s' , stb               , ...
                 'liqvol_r' , stb               , ...
                 'gasvol_s' , 1000 * ft^3       , ... % Mscf
                 'gasvol_r' , stb               , ...
                 'trans'    , centi*poise * stb / (day * psia));
      n = struct('length'   , 'ft'                , ...
                 'time'     , 'day'               , ...
                 'density'  , 'pound / ft^3'      , ...
                 'press'    , 'psia'              , ...
                 'concentr' , 'pound / stb'       , ... % Concentration
                 'compr'    , '1 / psia'          , ...
                 'viscosity', 'centi*poise'       , ...
                 'perm'     , 'milli*darcy'       , ...
                 'liqvol_s' , 'stb'               , ...
                 'liqvol_r' , 'stb'               , ...
                 'gasvol_s' , 'Mscf'       , ... % Mscf
                 'gasvol_r' , 'stb'               , ...
                 'trans'    , 'centi*poise * stb / (day * psia)');
   elseif lab,
      u = struct('length'   , centi*meter           , ...
                 'time'     , hour                  , ...
                 'density'  , gram / (centi*meter)^3, ...
                 'press'    , atm                   , ...
                 'concentr' , gram / (centi*meter)^3, ...
                 'compr'    , 1 / atm               , ...
                 'viscosity', centi*poise           , ...
                 'perm'     , milli*darcy           , ...
                 'liqvol_s' , (centi*meter)^3       , ...
                 'liqvol_r' , (centi*meter)^3       , ...
                 'gasvol_s' , (centi*meter)^3       , ...
                 'gasvol_r' , (centi*meter)^3       , ...
                 'trans'    , centi*poise * (centi*meter)^3 / (hour * atm));
      n = struct('length'   , 'centi*meter'           , ...
                 'time'     , 'hour'                  , ...
                 'density'  , 'gram / (centi*meter)^3', ...
                 'press'    , 'atm'                  , ...
                 'concentr' , 'gram / (centi*meter)^3', ...
                 'compr'    , '1 / atm'               , ...
                 'viscosity', 'centi*poise'           , ...
                 'perm'     , 'milli*darcy'           , ...
                 'liqvol_s' , '(centi*meter)^3'       , ...
                 'liqvol_r' , '(centi*meter)^3'       , ...
                 'gasvol_s' , '(centi*meter)^3'       , ...
                 'gasvol_r' , '(centi*meter)^3'       , ...
                 'trans'    , 'centi*poise * (centi*meter)^3 / (hour * atm)');
             
   else
      % SI units.  MRST extension.  Idempotency.
      u = struct('length'   , 1, ...
                 'time'     , 1, ...
                 'density'  , 1, ...
                 'press'    , 1, ...
                 'concentr' , 1, ...
                 'compr'    , 1, ...
                 'viscosity', 1, ...
                 'perm'     , 1, ...
                 'liqvol_s' , 1, ...
                 'liqvol_r' , 1, ...
                 'gasvol_s' , 1, ...
                 'gasvol_r' , 1, ...
                 'trans'    , 1);
      n = struct('length'   , 'meter', ...
                 'time'     , 'second', ...
                 'density'  , 'meter^3/kilogram', ...
                 'press'    , 'Pascal', ...
                 'concentr' , 'kilogram / meter^3', ...
                 'compr'    , '1 / Pascal'               , ...
                 'viscosity', 'Pascal*second'           , ...
                 'perm'     , 'meter^3'           , ...
                 'liqvol_s' , 'meter^3'       , ...
                 'liqvol_r' , 'meter^3'       , ...
                 'gasvol_s' , 'meter^3'       , ...
                 'gasvol_r' , 'meter^3'       , ...
                 'trans'    , 'meter^3');
   end
end