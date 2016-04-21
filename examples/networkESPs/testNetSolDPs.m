% REservoir Multiple Shooting Optimization.
% REduced Multiple Shooting Optimization.

% Make sure the workspace is clean before we start
    clc
    clear
    clear global

    % Required MRST modules
    mrstModule add deckformat
    mrstModule add ad-fi ad-core ad-props

    % Include REMSO functionalities
    addpath(genpath('../../mrstDerivated'));
    addpath(genpath('../../mrstLink'));
    addpath(genpath('../../mrstLink/wrappers/procedural'));
    
    addpath(genpath('../../netLink'));
    addpath(genpath('../../netLink/plottings'));
    addpath(genpath('../../netLink/fluidProperties'));
    addpath(genpath('../../netLink/networkFunctions'));
    addpath(genpath('../../netLink/auxiliaryFunctions'));

    addpath(genpath('../../optimization/multipleShooting'));
    addpath(genpath('../../optimization/plotUtils'));
    addpath(genpath('../../optimization/remso'));
    addpath(genpath('../../optimization/remsoSequential'));
    addpath(genpath('../../optimization/remsoCrossSequential'));
    addpath(genpath('../../optimization/singleShooting'));
    addpath(genpath('../../optimization/utils'));
    addpath(genpath('reservoirData'));


    %% Initialize reservoir -  the Simple reservoir
    [reservoirP] = initReservoir('RATE10x5x10.txt', 'Verbose',true);
%     [reservoirP] = initReservoir( 'reallySimpleRes.data','Verbose',true);

    % do not display reservoir simulation information!
    
    %% Greedy Algorithm Settings    
%     greedySteps = 5;    
%     controlSteps = arrayfun(@(x) numel(find(reservoirP.schedule.step.control==x)), 1:greedySteps);
%     
%     reservoirP.schedule.step.val =400*day;
%     reservoirP.schedule.step.control = 1;
%     reservoirP.schedule.control = reservoirP.schedule.control(1);     
%     
%     optYears = 4000*day;
%     greedySteps = optYears / reservoirP.schedule.step.val;        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
    mrstVerbose off;

    % Number of reservoir grid-blocks
    nCells = reservoirP.G.cells.num;

    %% Multiple shooting problem set up
    totalPredictionSteps = numel(reservoirP.schedule.step.val);  % MS intervals

    % Schedule partition for each control period and for each simulated step
    lastControlSteps = findControlFinalSteps( reservoirP.schedule.step.control );
    controlSchedules = multipleSchedules(reservoirP.schedule,lastControlSteps);

    stepSchedules = multipleSchedules(reservoirP.schedule,1:totalPredictionSteps);

    % Piecewise linear control -- mapping the step index to the corresponding
    % control
    ci  = arroba(@controlIncidence, 2 ,{reservoirP.schedule.step.control});    
    
    %% Variables Scaling
    xScale = setStateValues(struct('pressure',5*barsa,'sW',0.01),'nCells',nCells);

    if (isfield(reservoirP.schedule.control,'W'))
        W =  reservoirP.schedule.control.W;
    else
        W = processWells(reservoirP.G, reservoirP.rock,reservoirP.schedule.control(1),'DepthReorder', true);
    end
    wellSol = initWellSolLocal(W, reservoirP.state);
    for k = 1:numel(wellSol)
        wellSol(k).qGs = 0;
    end
    nW = numel(W);
    
    %% Fixing Injectors
%     fixedWells = find(vertcat(W.sign) == 1); % fixing injectors
    fixedWells = [];
    controlWells = setdiff(1:nW, fixedWells); 

    % Instantiate the production network object
    netSol = prodNetwork(wellSol, 'espNetwork', true);

    %%TODO: separate scalling of vk and nk.
    %% Scallings 
    [vScale, freqScale] = mrstAlg2algVar( wellSolScaling(wellSol,'bhp',5*barsa,'qWs',10*meter^3/day,'qOs',10*meter^3/day), netSolScaling(netSol));       
%     freqScale = [];
%     flowScale = [];        
    freqScale = [15;15;15;15;15]; % in Hz    
    flowScale = [5*(meter^3/day); 5*(meter^3/day);5*(meter^3/day);5*(meter^3/day);5*(meter^3/day); ...
                 5*(meter^3/day); 5*(meter^3/day);5*(meter^3/day);5*(meter^3/day);5*(meter^3/day)];

    pressureScale = [5*barsa;5*barsa;5*barsa;5*barsa;5*barsa];
       
%     pressureScale = [];
%     freqScale = [];
%     flowScale = [];
    
    %% network controls
    pScale = [];
    p  = [];
    
    % number of pump stages
    numStages =  [70; 80; 70; 80; 70];    
    % bounds for flowing rates through the pump at 60 Hz
    qlMin = [30*(meter^3/day); 30*(meter^3/day); 30*(meter^3/day); 30*(meter^3/day); 30*(meter^3/day)];
    qlMax = [150*(meter^3/day); 150*(meter^3/day); 150*(meter^3/day); 150*(meter^3/day); 150*(meter^3/day)];
    
    % bounds for pump frequencies
    freqMin = [30; 30; 30; 30; 30]; % in Hz
    freqMax = [90; 90; 90; 90; 90]; % in Hz
    baseFreq = [60; 60; 60; 60; 60]; % in Hz  
    
    
    %% run network four times to obtain the extreme points of the pumps maps        
    qminFmin = pump_rate(freqMin, qlMin, baseFreq);
    qminFmax = pump_rate(freqMax, qlMin, baseFreq);
    qmaxFmin = pump_rate(freqMin, qlMax, baseFreq);
    qmaxFmax = pump_rate(freqMax, qlMax, baseFreq);   
    
    qf = cell(4,1);
    dp = cell(4,1);
    
    qf{1}  = qminFmin./(meter^3/day); qf{2} = qminFmax./(meter^3/day); 
    qf{3} = qmaxFmin./(meter^3/day); qf{4} = qmaxFmax./(meter^3/day);            
    
    str = netSol.E(1).stream; % network has default stream for subsea pipeline
    
    oilDens = str.oil_dens;        
    highestDens = 0.6*str.water_dens + 0.4*str.oil_dens; % maximum of 0.6 water cut (based on experiments, not limited in practice)
    
    dp{1} = calcDp(qminFmin, freqMin, baseFreq, numStages, 'mixDensity',  oilDens);
    dp{2} = calcDp(qminFmax, freqMax, baseFreq, numStages, 'mixDensity', highestDens);
    dp{3} = calcDp(qmaxFmin, freqMin, baseFreq, numStages, 'mixDensity', oilDens);
    dp{4} = calcDp(qmaxFmax, freqMax, baseFreq, numStages, 'mixDensity', highestDens);
    
    extremePoints = cell(4,1);
    for i=1:4
        extremePoints{i} = [qf{i}, dp{i}];        
    end
    extremePoints = [];

    
    % function that performs a network simulation, and calculates the
    % pressure drop (dp) in the chokes/pumps  
    dpPumps = arroba(@pumpsDp,[1,2,3],{netSol, pressureScale, numStages, pScale, 'turnoffPumps', false}, true);        
    pumpFrequencies = arroba(@pumpFrequency,[1,2,3],{netSol, freqScale, numStages, pScale}, true);    
    minFlowPump = arroba(@pumpFlowMin,[1,2,3],{netSol, flowScale, numStages, pScale}, true);
    maxFlowPump = arroba(@pumpFlowMax,[1,2,3],{netSol, flowScale, numStages, pScale}, true);
      

    cellControlScales = schedules2CellControls(schedulesScaling(controlSchedules,...
        'RATE',10*meter^3/day,...
        'ORAT',10*meter^3/day,...
        'WRAT',10*meter^3/day,...
        'LRAT',10*meter^3/day,...
        'RESV',0,...
        'BHP',5*barsa),'fixedWells', fixedWells);

    %% instantiate the objective function as an aditional Algebraic variable


    %%% The sum of the last elements in the algebraic variables is the objective
    nCells = reservoirP.G.cells.num;
    stepNPV = arroba(@NPVStepM,[1,2, 3],{nCells,'scale',1/100000,'sign',-1},true);    

    nScale  = [flowScale; freqScale; pressureScale];  
   %nScalePump = [freqScale; pressureScale];  
   
    vScale = [vScale; nScale; 1];
	%vScalePump = [vScale; nScalePump; 1];
    
    networkJointObj = arroba(@networkJointNPVConstraints,[1,2, 3],{nCells, netSol, freqScale, pressureScale, flowScale, numStages, baseFreq, qlMin, qlMax, pScale,   'scale',1/100000,'sign',-1, 'turnoffPumps', false, 'dpFunction', @simpleDp, 'extremePoints', extremePoints},true);

%     [ algFun ] = concatenateMrstTargets(networkJointObj,false, [numel(nScale); 1]);
    %[ algFunPump ] = concatenateMrstTargets([pumpFrequencies, stepNPV],false, [numel(nScalePump); 1]);
    
%     [ algFun ] = concatenateMrstTargets([stepNPV],false, [numel(vScale); 1]);
    
%     [ algFun ] = concatenateMrstTargets([pumpFrequencies, dpPumps, maxFlowPump, stepNPV],false, [numel(vScale); 1]);

%     [ algFun ] = concatenateMrstTargets([pumpFrequencies, dpPumps, stepNPV],false, [numel(vScale); 1]);

    %% Instantiate the simulators for each interval, locally and for each worker.

    % ss.stepClient = local (client) simulator instances
    % ss.state = scaled initial state
    % ss.nv = number of algebraic variables
    % ss.ci = piecewice control mapping on the client side
    % ss.step =  worker simulator instances
    
    stepBreak = 40;

    step = cell(totalPredictionSteps,1);                                            
    for k=1:totalPredictionSteps
        cik = callArroba(ci,{k});
%         if k<=stepBreak
%             step{k} = @(x0,u,varargin) mrstStep(x0,u,@mrstSimulationStep,wellSol,stepSchedules(k),reservoirP,...
%                                             'xScale',xScale,...
%                                             'vScale',vScale,...
%                                             'uScale',cellControlScales{cik},...
%                                             'algFun',algFunPump,...
%                                             'fixedWells', fixedWells, ...
%                                             varargin{:});
%         else
            step{k} = @(x0,u,varargin) mrstStep(x0,u,@mrstSimulationStep,wellSol,stepSchedules(k),reservoirP,...
                                            'xScale',xScale,...
                                            'vScale',vScale,...
                                            'uScale',cellControlScales{cik},...
                                            'algFun',networkJointObj,...
                                            'fixedWells', fixedWells, ...
                                            varargin{:});
%         end
                                        
    end

    ss.state = stateMrst2stateVector( reservoirP.state,'xScale',xScale );
    ss.step = step;
    ss.ci = ci;

    %% instantiate the objective function
    %%% objective function
    obj = cell(totalPredictionSteps,1);
    for k = 1:totalPredictionSteps
        obj{k} = arroba(@lastAlg,[1,2,3],{},true);
    end
    targetObj = @(xs,u,vs,varargin) sepTarget(xs,u,vs,obj,ss,varargin{:});

    %%  Bounds for all variables!

    % Bounds for all wells!
    % minProd = struct('BHP',130*barsa, 'ORAT', 1*meter^3/day); original val
    minProd = struct('BHP',100*barsa, 'ORAT',  5*meter^3/day);

    % maxProd = struct('BHP',200*barsa, 'ORAT', 220*meter^3/day); original val
    maxProd = struct('BHP',500*barsa, 'ORAT', 200*meter^3/day); 

    % minInj = struct('RATE',100*meter^3/day); % original val
    minInj = struct('RATE',1*meter^3/day);
    % maxInj = struct('RATE',300*meter^3/day); original val

    maxInj = struct('RATE',300*meter^3/day);

    % Control input bounds for all wells!

    [ lbSchedules,ubSchedules ] = scheduleBounds( controlSchedules,...
        'maxProd',maxProd,'minProd',minProd,...
        'maxInj',maxInj,'minInj',minInj,'useScheduleLims',false);
    lbw = schedules2CellControls(lbSchedules,'cellControlScales',cellControlScales, 'fixedWells', fixedWells);
    ubw = schedules2CellControls(ubSchedules,'cellControlScales',cellControlScales, 'fixedWells', fixedWells);   


    cellControlScale = cellfun(@(wi) [wi; pScale],cellControlScales,'uniformOutput', false);

    lbu = cellfun(@(wi)[wi; 5*barsa./pScale],lbw, 'UniformOutput',false);
    ubu = cellfun(@(wi)[wi; 30*barsa./pScale],ubw, 'UniformOutput',false);


    % Bounds for all wells!
    % minProd = struct('ORAT',1*meter^3/day,  'WRAT',1*meter^3/day,  'GRAT',
    % -inf,'BHP',130*barsa); original val
    minProd = struct('ORAT', 5*meter^3/day,  'WRAT', 0*meter^3/day,  'GRAT', -inf,'BHP',100*barsa);
    % maxProd = struct('ORAT',220*meter^3/day,'WRAT',150*meter^3/day,'GRAT',
    % inf,'BHP',350*barsa); original val
    maxProd = struct('ORAT',200*meter^3/day,'WRAT',200*meter^3/day,'GRAT', inf,'BHP',500*barsa);

    % minInj = struct('ORAT',-inf,  'WRAT',100*meter^3/day,  'GRAT',
    % -inf,'BHP', 5*barsa); original val
    minInj = struct('ORAT',-inf,  'WRAT', 1*meter^3/day,  'GRAT', -inf,'BHP', 5*barsa);
    % maxInj = struct('ORAT',inf,'WRAT',300*meter^3/day,'GRAT',
    % inf,'BHP',500*barsa); original val
    maxInj = struct('ORAT',inf,'WRAT', 300*meter^3/day,'GRAT', inf,'BHP',800*barsa);

    % wellSol bounds  (Algebraic variables bounds)
    [ubWellSol,lbWellSol] = wellSolScheduleBounds(wellSol,...
        'maxProd',maxProd,...
        'maxInj',maxInj,...
        'minProd',minProd,...
        'minInj',minInj);

    ubvS = wellSol2algVar(ubWellSol,'vScale',vScale);
    lbvS = wellSol2algVar(lbWellSol,'vScale',vScale);    

    %% Linear Approx. of Pump Map
%     lbv = repmat({[lbvS; 0./flowScale;   0*barsa./pressureScale; 0*barsa./pressureScale;  -inf*barsa./pressureScale; -inf]},totalPredictionSteps,1);
%     ubv = repmat({[ubvS; inf./flowScale; inf*barsa./pressureScale; inf*barsa./pressureScale;  inf*barsa./pressureScale; inf]},totalPredictionSteps,1);

    %% Non-Linear Pump Map Constraints
    lbv = repmat({[lbvS; 0./flowScale;   freqMin./freqScale;  -inf*barsa./pressureScale; -inf]},totalPredictionSteps,1);
    ubv = repmat({[ubvS; inf./flowScale; freqMax./freqScale;  inf*barsa./pressureScale; inf]},totalPredictionSteps,1);


   %% Unconstrained problem regarding the network
%     lbv = repmat({[lbvS; -inf]},totalPredictionSteps,1);
%     ubv = repmat({[ubvS;  inf]},totalPredictionSteps,1);      

    % State lower and upper - bounds
    maxState = struct('pressure',800*barsa,'sW',1);
    minState = struct('pressure',50*barsa,'sW',0.1);
    ubxS = setStateValues(maxState,'nCells',nCells,'xScale',xScale);
    lbxS = setStateValues(minState,'nCells',nCells,'xScale',xScale);
    lbx = repmat({lbxS},totalPredictionSteps,1);
    ubx = repmat({ubxS},totalPredictionSteps,1);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    lowActive = [];
    upActive = [];

    %% A plot function to display information at each iteration

    times.steps = [stepSchedules(1).time;arrayfun(@(x)(x.time+sum(x.step.val))/day,stepSchedules)];
    times.tPieceSteps = cell2mat(arrayfun(@(x)[x;x],times.steps,'UniformOutput',false));
    times.tPieceSteps = times.tPieceSteps(2:end-1);

    times.controls = [controlSchedules(1).time;arrayfun(@(x)(x.time+sum(x.step.val))/day,controlSchedules)];
    times.tPieceControls = cell2mat(arrayfun(@(x)[x;x],times.controls,'UniformOutput',false));
    times.tPieceControls = times.tPieceControls(2:end-1);

    cellControlScalesPlot = schedules2CellControls(schedulesScaling( controlSchedules,'RATE',1/(meter^3/day),...
        'ORAT',1/(meter^3/day),...
        'WRAT',1/(meter^3/day),...
        'LRAT',1/(meter^3/day),...
        'RESV',0,...
        'BHP',1/barsa));

    cellControlScalesPlot = cellfun(@(w) [w;pScale], cellControlScalesPlot, 'UniformOutput',false); 

    cellControlScales  = cellfun(@(w) [w; pScale] , cellControlScales ,'uniformOutput', false);

    [uMlb] = scaleSchedulePlot(lbu,controlSchedules,cellControlScales,cellControlScalesPlot, 'fixedWells', fixedWells);
    [uLimLb] = min(uMlb,[],2);
    ulbPlob = cell2mat(arrayfun(@(x)[x,x],uMlb,'UniformOutput',false));


    [uMub] = scaleSchedulePlot(ubu,controlSchedules,cellControlScales,cellControlScalesPlot, 'fixedWells', fixedWells);
    [uLimUb] = max(uMub,[],2);
    uubPlot = cell2mat(arrayfun(@(x)[x,x],uMub,'UniformOutput',false));


    % be carefull, plotting the result of a forward simulation at each
    % iteration may be very expensive!
    % use simFlag to do it when you need it!
    simFunc =@(sch,varargin) runScheduleADI(reservoirP.state, reservoirP.G, reservoirP.rock, reservoirP.system, sch,'force_step',false,varargin{:});

    wc    = vertcat(W.cells);
    fPlot = @(x)[max(x);min(x);x(wc)];

    %prodInx  = (vertcat(wellSol.sign) < 0);
    %wc    = vertcat(W(prodInx).cells);
    %fPlot = @(x)x(wc);

    % plotSol = @(x,u,v,d,varargin) plotSolution( x,u,v,d,ss,obj,times,xScale,cellControlScales,vScale,cellControlScalesPlot,controlSchedules,wellSol,ulbPlob,uubPlot,[uLimLb,uLimUb],minState,maxState,'simulate',simFunc,'plotWellSols',true,'plotSchedules',false,'pF',fPlot,'sF',fPlot,varargin{:});

    plotSol = @(x,u,v,d,varargin) plotSolution( x,u,v,d, lbv, ubv, lbu, ubu, ss,obj,times,xScale,cellControlScales,vScale, nScale, ...
                cellControlScalesPlot,controlSchedules,wellSol,ulbPlob,uubPlot,[uLimLb,uLimUb],minState,maxState,'simulate',simFunc,'plotWellSols',true, 'plotNetsol', true, ...
                'numNetConstraints', numel(nScale), 'plotNetControls', false, 'numNetControls', numel(pScale), 'freqCst', numel(freqScale), 'pressureCst',numel(pressureScale),  'flowCst',numel(flowScale), ...
                'plotSchedules',false,'pF',fPlot,'sF',fPlot, 'fixedWells', fixedWells, 'extremePoints', extremePoints, 'plotCumulativeObjective', true, 'qlMin', qlMin,  'qlMax', qlMax, 'nStages', numStages, ...
                'freqMin', freqMin, 'freqMax', freqMax, 'baseFreq', baseFreq, varargin{:});
            

    % remove network control to initialize well controls vector (w)
    cellControlScales = cellfun(@(w) w(1:end-numel(p)) ,cellControlScales, 'UniformOutput', false);

    %%  Initialize from previous solution?

    x = [];
    v = [];
    w  = schedules2CellControls( controlSchedules,'cellControlScales',cellControlScales, 'fixedWells', fixedWells);
   
    cellControlScales = cellfun(@(w) [w; pScale] , cellControlScales ,'uniformOutput', false);
    u = cellfun(@(wi)[wi;p],w,'UniformOutput',false);

    cellControlScalesPlot = cellfun(@(w) [w; pScale], cellControlScalesPlot,'uniformOutput', false);

    testFlag = false;
    if testFlag    
        addpath(genpath('../../optimization/testFunctions'));

        [~, ~, ~, simVars, xs, vs] = simulateSystemSS(u, ss, []);
        [ei, fi, vi] = testProfileGradients(xs,u,vs,ss.step,ss.ci,ss.state, 'd', 1, 'pert', 1e-5, 'all', false);       

    end
    
    
    %%%% Unit test for pressure drop calculation
    %%%% 

    warning('ignoring gas')
    
    nProd = sum(vertcat(W.sign)==-1);
    qoBounds = [minProd.ORAT;maxProd.ORAT*nProd];
    qwBounds = [minProd.WRAT;maxProd.WRAT*nProd];
    
    warning('ignoring gas')
    qgBounds = [0;0];  % there is no gas in the problem
 
    sepVertex = getVertex(netSol, setdiff(netSol.Vsnk, netSol.VwInj));
    pBounds = [sepVertex.pressure;maxProd.BHP];    
    
    nGrid = 100;
    
    qoGrid = qoBounds(1):(qoBounds(2)-qoBounds(1))/nGrid:qoBounds(2);
    qwGrid = qwBounds(1):(qwBounds(2)-qwBounds(1))/nGrid:qwBounds(2);
    qgGrid = qgBounds(1):(qgBounds(2)-qgBounds(1))/nGrid:qgBounds(2);
    pGrid  =  pBounds(1):( pBounds(2)- pBounds(1))/nGrid: pBounds(2);
    
    pert = 1e-5;
    qoP = 5*meter^3/day*pert;
    qwP = 5*meter^3/day*pert;
    qgP = 100*(10*ft)^3/day*pert;
    pP =  5*barsa*pert;
    
    gradScale = 5*barsa./[5*meter^3/day,5*meter^3/day,100*(10*ft)^3/day,5*barsa];

    if isempty(qoGrid)
        qoGrid = qoBounds(1);
    end    
    if isempty(qwGrid)
        qwGrid = qwBounds(1);
    end    
    if isempty(qgGrid)
        qgGrid = qgBounds(1);
    end
    if isempty(pGrid)
        pGrid = pBounds(1);
    end    
    
    wcutPlot = true;
    testGradient = true;
    % iterator for all pipelines
    for i=1:numel(netSol.E)
        ei = netSol.E(i);
        
        
        if wcutPlot
            %figure(1);hold on;
            
            %qlGrid = (logspace(log(0.001)/log(10),log(200)/log(10),100))*meter^3/day;
            qlMax = 200;
            n = 100;
            qlGrid = (0.001:(qlMax-0.001)/n:qlMax)*meter^3/day;
            wcutGrid = 0:0.25:1;
            qgGrid = 0;
            pGrid = 150*barsa;  % dpBeggsBrill is insensitive to pV when there is no gas;

             %(1)  simulate for the grid and plot
            dpGrid = zeros(numel(qlGrid),numel(wcutGrid ),numel(qgGrid),numel(pGrid));
            for pk = 1:numel(pGrid)   
                for gk = 1:numel(qgGrid)
                    for wk = 1:numel(wcutGrid )                        
                        for ok = 1:numel(qlGrid) 
                            dpGrid(ok,wk,gk,pk) = dpCVODES(ei,  qlGrid(ok)*(1-wcutGrid(wk)), qlGrid(ok)*wcutGrid(wk), qgGrid(gk), pGrid(pk), 'dpFunction', @dpBeggsBrill);
                        end
            %            plot(qlGrid*day,dpGrid(:,wk,gk,pk)/barsa,'.-')
                    end 
                end
            end
            close all;
            figure(); hold on;grid
            for wk = 1:numel(wcutGrid )
                plot(qlGrid/(meter^3/day),dpGrid(:,wk,gk,pk)/barsa,'.-','Color',[0,0,wcutGrid(wk)])
            end
            title(ei.name)
            ylabel('dp (bar)')
            xlabel('Liquid rate (meter^3/day)')
            annotation('textbox',[0.12,0.92,0,0],...
            'String',{['diam=' num2str(ei.pipeline.diam),'m.'],...
                      ['leng=' num2str(ei.pipeline.len),'m.'],...
                      ['angl=' num2str(rad2deg(ei.pipeline.ang)),'deg.'],...
                      ['temp=' num2str(ei.pipeline.temp),'K.']});
            savefig([ei.name,'_wCutGrid']);
            saveas(gcf,[ei.name,'_wCutGrid','.png']);
            save([ei.name,'_wCutGrid'],'ei','qlGrid','wcutGrid','qgGrid','pGrid','dpGrid');       
        
        else

            %(1)  simulate for the grid and plot
            dpGrid = zeros(numel(qoGrid),numel(qwGrid),numel(qgGrid),numel(pGrid));
            for ok = 1:numel(qoGrid)
                for wk = 1:numel(qwGrid)
                    for gk = 1:numel(qgGrid)
                        for pk = 1:numel(pGrid)                        
                            dpGrid(ok,wk,gk,pk) = dpCVODES(ei,  qoGrid(ok), qwGrid(wk), qgGrid(gk), pGrid(pk), 'dpFunction', @dpBeggsBrill);
                        end
                    end
                end
            end

            save([ei.name,'_grid'],qoGrid,qwGrid,qgGrid,pGrid,dpGrid);
        
        end
        
        % save the plot data to see afterwards.
        
        %(2) Pick random numbers within the bounds test gradients againts FD
        if testGradient

            nGrid = 1;
            qoRandom = qoBounds(1) + (qoBounds(2)-qoBounds(1))*rand(nGrid,1);
            qwRandom = qwBounds(1) + (qwBounds(2)-qwBounds(1))*rand(nGrid,1);
            qgRandom = qgBounds(1) + (qgBounds(2)-qgBounds(1))*rand(nGrid,1);
             pRandom =  pBounds(1) + ( pBounds(2)- pBounds(1))*rand(nGrid,1);



            dpErrAbs = zeros(nGrid,1);
            for k = 1:nGrid
                [qok, qwk, qgk,pk] = deal(qoRandom(k), qwRandom(k), qgRandom(k), pRandom(k));

                dpGradNum =[ 
                            (dpCVODES(ei, qok+qoP, qwk, qgk,pk , 'dpFunction', @dpBeggsBrill)-dpCVODES(ei, qok-qoP, qwk, qgk,pk , 'dpFunction', @dpBeggsBrill))/(2*qoP),...
                            (dpCVODES(ei, qok, qwk+qwP, qgk,pk , 'dpFunction', @dpBeggsBrill)-dpCVODES(ei, qok, qwk-qwP, qgk,pk , 'dpFunction', @dpBeggsBrill))/(2*qwP),...
                            0,...
                            (dpCVODES(ei, qok, qwk, qgk,pk+ pP , 'dpFunction', @dpBeggsBrill)-dpCVODES(ei, qok, qwk, qgk,pk- pP , 'dpFunction', @dpBeggsBrill))/(2* pP),...
                           ];


                [qok, qwk, qgk,pk] = initVariablesADI(qoRandom(k), qwRandom(k), qgRandom(k), pRandom(k)); 
                dpk = dpCVODES(ei, qok, qwk, qgk,pk , 'dpFunction', @dpBeggsBrill);
                dpGrad = cell2mat(dpk.jac);

                dpErrAbs(k) =norm(   (dpGrad-dpGradNum)./gradScale  )
            end       

            save([ei.name,'gradRand'],'qoRandom','qwRandom','qgRandom','pRandom','dpErrAbs');
        
        end
        
    end
    
    
    
    