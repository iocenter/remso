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
    
    addpath(genpath('../../optimization/multipleShooting'));
    addpath(genpath('../../optimization/plotUtils'));
    addpath(genpath('../../optimization/remso'));
    addpath(genpath('../../optimization/remsoSequential'));
    addpath(genpath('../../optimization/remsoCrossSequential'));
    addpath(genpath('../../optimization/singleShooting'));
    addpath(genpath('../../optimization/utils'));
    addpath(genpath('reservoirData'));


    %% Initialize reservoir -  the Simple reservoir
    [reservoirP] = initReservoir('SIMPLE10x5x10.txt', 'Verbose',true);
%     [reservoirP] = initReservoir( 'reallySimpleRes.data','Verbose',true);

    % do not display reservoir simulation information!
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
    fixedWells = find(vertcat(W.sign) == 1); % fixing injectors
    controlWells = setdiff(1:nW, fixedWells); 

    % Instantiate the production network object
    netSol = prodNetwork(wellSol, 'espNetwork', true);

    %%TODO: separate scalling of vk and nk.
    %% Scallings 
    [vScale, freqScale] = mrstAlg2algVar( wellSolScaling(wellSol,'bhp',5*barsa,'qWs',10*meter^3/day,'qOs',10*meter^3/day), netSolScaling(netSol));       
%     freqScale = [];
%     flowScale = [];        
%     freqScale = [15;15;15;15;15]; % in Hz    
%     flowScale = [5*(meter^3/day); 5*(meter^3/day);5*(meter^3/day);5*(meter^3/day);5*(meter^3/day)];
    freqScale = [];
    flowScale = [];
    pressureScale = [5*barsa;5*barsa;5*barsa;5*barsa;5*barsa];
       
    
    %% network controls
    pScale = [];
    p  = [];
    
    % number of pump stages
    numStages =  30./ones(numel(freqScale),1);

    % function that performs a network simulation, and calculates the
    % pressure drop (dp) in the chokes/pumps  
    %% TODO: correct jacobian size corresponding to fixed wells
  
    dpPumps = arroba(@pumpsDp,[1,2,3],{netSol, pressureScale, numStages, pScale, 'turnoffPumps', true}, true);        
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
    stepGreedy = arroba(@greedyStep,[1,2, 3],{nCells,'scale',1/100000,'sign',-1},true);

    nScale  = [freqScale; pressureScale; flowScale];  
   
    vScale = [vScale; nScale; 1];

    [ algFun ] = concatenateMrstTargets([dpPumps, stepGreedy],false, [numel(nScale); 1]);
    
%     [ algFun ] = concatenateMrstTargets([dpPumps, stepNPV],false, [numel(vScale); 1]);
    
%     [ algFun ] = concatenateMrstTargets([pumpFrequencies, dpPumps, maxFlowPump, stepNPV],false, [numel(vScale); 1]);

%     [ algFun ] = concatenateMrstTargets([pumpFrequencies, dpPumps, stepNPV],false, [numel(vScale); 1]);

    %% Instantiate the simulators for each interval, locally and for each worker.

    % ss.stepClient = local (client) simulator instances
    % ss.state = scaled initial state
    % ss.nv = number of algebraic variables
    % ss.ci = piecewice control mapping on the client side
    % ss.step =  worker simulator instances
    


    step = cell(totalPredictionSteps,1);
    for k=1:totalPredictionSteps
        cik = callArroba(ci,{k});
        step{k} = @(x0,u,varargin) mrstStep(x0,u,@mrstSimulationStep,wellSol,stepSchedules(k),reservoirP,...
                                            'xScale',xScale,...
                                            'vScale',vScale,...
                                            'uScale',cellControlScales{cik},...
                                            'algFun',algFun,...
                                            'fixedWells', fixedWells, ...
                                            varargin{:});
                                        
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
    minProd = struct('BHP',1*barsa, 'ORAT',  5*meter^3/day);

    % maxProd = struct('BHP',200*barsa, 'ORAT', 220*meter^3/day); original val
    maxProd = struct('BHP',800*barsa, 'ORAT', 500*meter^3/day); 

    % minInj = struct('RATE',100*meter^3/day); % original val
    minInj = struct('RATE',1*meter^3/day);
    % maxInj = struct('RATE',300*meter^3/day); original val

    maxInj = struct('RATE',500*meter^3/day);

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
    minProd = struct('ORAT', 5*meter^3/day,  'WRAT', -inf*meter^3/day,  'GRAT', -inf,'BHP',5*barsa);
    % maxProd = struct('ORAT',220*meter^3/day,'WRAT',150*meter^3/day,'GRAT',
    % inf,'BHP',350*barsa); original val
    maxProd = struct('ORAT',500*meter^3/day,'WRAT',500*meter^3/day,'GRAT', inf,'BHP',800*barsa);

    % minInj = struct('ORAT',-inf,  'WRAT',100*meter^3/day,  'GRAT',
    % -inf,'BHP', 5*barsa); original val
    minInj = struct('ORAT',-inf,  'WRAT', 1*meter^3/day,  'GRAT', -inf,'BHP', 5*barsa);
    % maxInj = struct('ORAT',inf,'WRAT',300*meter^3/day,'GRAT',
    % inf,'BHP',500*barsa); original val
    maxInj = struct('ORAT',inf,'WRAT', 500*meter^3/day,'GRAT', inf,'BHP',800*barsa);

    % wellSol bounds  (Algebraic variables bounds)
    [ubWellSol,lbWellSol] = wellSolScheduleBounds(wellSol,...
        'maxProd',maxProd,...
        'maxInj',maxInj,...
        'minProd',minProd,...
        'minInj',minInj);

    ubvS = wellSol2algVar(ubWellSol,'vScale',vScale);
    lbvS = wellSol2algVar(lbWellSol,'vScale',vScale);

    nScale  = [freqScale; pressureScale; flowScale];
    
    lbv = repmat({[lbvS; -inf*barsa./pressureScale; -inf]},totalPredictionSteps,1);
    ubv = repmat({[ubvS;  inf*barsa./pressureScale; inf]},totalPredictionSteps,1);
    
    % constrain pump by the pressure loss given by the max flow rate.
%      lbv = repmat({[lbvS; 20 ./freqScale;  -100*barsa./pressureScale;-inf]},totalPredictionSteps,1);
%      ubv = repmat({[ubvS; 100./freqScale; 100*barsa./pressureScale;inf]},totalPredictionSteps,1);
    
%     lbv = repmat({[lbvS; -inf./flowScale;  -inf./flowScale;  30./freqScale;  -inf]},totalPredictionSteps,1);
%     ubv = repmat({[ubvS;  0./flowScale; 0./flowScale; 60./freqScale;  inf]},totalPredictionSteps,1);

%     lbv = repmat({[lbvS; -inf./flowScale;  -inf./flowScale;  0./freqScale;  -inf]},totalPredictionSteps,1);
%     ubv = repmat({[ubvS;  inf./flowScale; inf./flowScale; inf./freqScale;  inf]},totalPredictionSteps,1);

%     lbv = repmat({[lbvS;  20./freqScale; -100*barsa./pressureScale; -inf./flowScale;  -inf]},totalPredictionSteps,1);
%     ubv = repmat({[ubvS;  100./freqScale; 100*barsa./pressureScale; 0./flowScale;  inf]},totalPredictionSteps,1);

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
                'plotSchedules',false,'pF',fPlot,'sF',fPlot, 'fixedWells', fixedWells, varargin{:});

    % remove network control to initialize well controls vector (w)
    cellControlScales = cellfun(@(w) w(1:end-numel(p)) ,cellControlScales, 'UniformOutput', false);

    %%  Initialize from previous solution?

    % if exist('optimalVars.mat','file') == 2
    %     load('optimalVars.mat','x','u','v');
    % elseif exist('itVars.mat','file') == 2
    %     load('itVars.mat','x','u','v');
    % else
        x = [];
        v = [];
        w  = schedules2CellControls( controlSchedules,'cellControlScales',cellControlScales, 'fixedWells', fixedWells);
        %[x] = repmat({ss.state},totalPredictionSteps,1);
    % end

    cellControlScales = cellfun(@(w) [w; pScale] , cellControlScales ,'uniformOutput', false);
    u = cellfun(@(wi)[wi;p],w,'UniformOutput',false);

    cellControlScalesPlot = cellfun(@(w) [w; pScale], cellControlScalesPlot,'uniformOutput', false);

    testFlag = false;
    if testFlag    
        addpath(genpath('../../optimization/testFunctions'));

        [~, ~, ~, simVars, xs, vs] = simulateSystemSS(u, ss, []);
        [ei, fi, vi] = testProfileGradients(xs,u,vs,ss.step,ss.ci,ss.state, 'd', 1, 'pert', 1e-5, 'all', true);
    end

algorithm = 'remso';
switch algorithm
    
    case 'remso'
        %% call REMSO
        
%         load itVars;       

        [u,x,v,f,xd,M,simVars] = remso(u,ss,targetObj,'lbx',lbx,'ubx',ubx,'lbv',lbv,'ubv',ubv,'lbu',lbu,'ubu',ubu,...
            'tol',1e-6,'lkMax',4,'debugLS',true,...
            'lowActive',lowActive,'upActive',upActive,...
            'plotFunc',plotSol,'max_iter',500,'x',x,'v',v,'debugLS',false,'saveIt',true, 'computeCrossTerm', false, 'condense', true);
        
        %% plotSolution     
%{        
         [~, ~, ~, simVars, x, v] = simulateSystemSS(u, ss, []);
         xd = cellfun(@(xi)xi*0,x,'UniformOutput',false);
         plotSol(x,u,v,xd, 'simFlag', false);        
%}  
         plotSol(x,u,v,xd, 'simFlag', false);    
         
    case 'greedy'
     
        %% loop controlling the interval (timestep)
        steps = reservoirP.schedule.step.val;
        controlSteps = reservoirP.schedule.step.control;
        controls     =  schedule.control;
        for i=1:numel(steps)
            
            reservoirP.schedule.step.val = steps(i);
            reservoirP.schedule.step.control = controlSteps(i);
            reservoirP.schedule.control = controls(controlSteps(i));
            
            %% call remso to optimize the interval
            [u,x,v,f,xd,M,simVars] = remso(u,ss,targetObj,'lbx',lbx,'ubx',ubx,'lbv',lbv,'ubv',ubv,'lbu',lbu,'ubu',ubu,...
                'tol',1e-6,'lkMax',4,'debugLS',true,...
                'lowActive',lowActive,'upActive',upActive,...
                'plotFunc',plotSol,'max_iter',500,'x',x,'v',v,'debugLS',false,'saveIt',true, 'computeCrossTerm', false, 'condense', true);
            
            
            
            save (i, 'p', 'q'); 
            
            %% save controls (solution) of the step (itVars.mat)
        end
        
        %% apply the controls in a forward simulation to obtain the results with the greedy strategy
         
    case 'snopt'
        
        objSparsity = ones(1,size(cell2mat(u),1));
        
        uDim = cellfun(@(x)size(x,1),u);
        [outputCons,lbC,ubC,consSparsity] = outputVarsBoundSelector(lbx,ubx,lbv,ubv,uDim,ci);
        
        consSizes = cellfun(@(x)size(x,1),lbC);
        cons = cell(numel(consSizes),1);
        for k=1:size(cons)
            cons{k} = arroba(@concatenateTargetK,[2,3,4],{k,outputCons{k},consSizes},true);
        end
        
        outDims = [1,sum(cellfun(@(x)size(x,1),consSparsity))];
        [ target ] = concatenateTargets(obj,cons,outDims);       
        
        
        objCons = @(u,varargin) simulateSystemSS(u,ss,target,'abortNotConvergent',true,varargin{:});
        
        
        objGradFG = @(uu) dealSnoptSimulateSS( uu,objCons,cellfun(@(x)numel(x),u),true);
        
        optionsSNOPT = which('options.spc');
        if ~strcmp(optionsSNOPT,'')
            snspec(optionsSNOPT);
        end
        
        sparsity = [objSparsity;cell2mat(consSparsity)];
        
        
        snset  ('Minimize');
        snseti('Derivative Option',1);
        snscreen on
        
        ObjAdd = 0;
        ObjRow = 1;
        
        A= [];
        iAfun = [];
        jAvar = [];
        
        [iGfun,jGvar] = find(sparsity);
        
        if size(iGfun,1) < size(iGfun,2)
            iGfun = iGfun';
            jGvar = jGvar';
        end
        
        if exist('snoptLog.txt','file')==2
            delete('snoptLog.txt');
        end
        if exist('snoptSummary.txt','file')==2
            delete('snoptSummary.txt');
        end
        if exist('snoptDetail.txt','file')==2
            delete('snoptDetail.txt');
        end
        snsummary( 'snoptSummary.txt');
        snprintfile( 'snoptDetail.txt');
        
        [u,F,inform,xmul,Fmul] = snopt(cell2mat(u),...
            cell2mat(lbu),...
            cell2mat(ubu),...
            [-inf;cell2mat(lbC)],...
            [ inf;cell2mat(ubC)],...
            objGradFG,...
            ObjAdd,ObjRow,...
            A, iAfun, jAvar, iGfun, jGvar);
        
        snsummary( 'off');
        snprintfile( 'off');
        
        %{
        u = u.*cell2mat(cellControlScales);
        uC = mat2cell(u,uDims,1);
        schedule = cellControls2Schedule( uC,reservoirP.schedule );
        [wellSols,States] = simFunc(schedule);
        [qWs, qOs, qGs, bhp] = wellSolToVector(wellSols);
        
        figure(1)
        plot(cumsum(reservoirP.schedule.step.val),qWs*day)
        title('water (meter^3/day)')
        
        figure(2)
        plot(cumsum(reservoirP.schedule.step.val),qOs*day)
        title('oil (meter^3/day)')
        
        figure(3)
        plot(cumsum(reservoirP.schedule.step.val),bhp/barsa)
        title('bhp (barsa)')
        %}
    case 'ipopt'
        objectiveSS = @(u,varargin) simulateSystemSS(u,ss,obj,varargin{:});
        
        objSparsity = ones(1,size(cell2mat(u),1));
        
        uDim = cellfun(@(x)size(x,1),u);
        [outputCons,lbC,ubC,consSparsity] = outputVarsBoundSelector(lbx,ubx,lbv,ubv,uDim,ci);
        
        consSizes = cellfun(@(x)size(x,1),lbC);
        cons = cell(numel(consSizes),1);
        for k=1:size(cons)
            cons{k} = @(xsk,vsk,uk,varargin) concatenateTargetK(k,xsk,vsk,uk,outputCons{k},consSizes,varargin{:});
        end
        
        constraintSS = @(u,varargin) simulateSystemSS(u,ss,cons,varargin{:});
        
        
        x0         = cell2mat(u);   % The starting point.
        options.lb = cell2mat(lbu);  % Lower bound on the variables.
        options.ub = cell2mat(ubu);  % Upper bound on the variables.
        options.cl = cell2mat(lbC);   % Lower bounds on the constraint functions.
        options.cu = cell2mat(ubC);   % Upper bounds on the constraint functions.
        
        [ fM ] = memorizeLastSimulation(u,[],true);
        
        
        fwdObj = @(x) ssFwdMemory(mat2cell(x,uDim,1),...
            @(xx,varargin)objectiveSS(xx,'gradients',false,varargin{:}),...
            fM,...
            'replace',false);
        gradObj = @(x) cell2mat(ssFwdMemory(mat2cell(x,uDim,1),...
            @(xx,varargin)objectiveSS(xx,'gradients',true,varargin{:}),...
            fM,...
            'replace',true,'gradFirst',true));
        fwdCons = @(x) ssFwdMemory(mat2cell(x,uDim,1),...
            @(xx,varargin)constraintSS(xx,'gradients',false,varargin{:}),...
            fM,...
            'replace',false);
        gradCons = @(x) sparse(cell2mat(ssFwdMemory(mat2cell(x,uDim,1),...
            @(xx,varargin)constraintSS(xx,'gradients',true,varargin{:}),...
            fM,...
            'replace',true,'gradFirst',true)));
        
        % The callback functions.
        funcs.objective        = fwdObj;
        funcs.gradient         = gradObj;
        funcs.constraints       = fwdCons;
        funcs.jacobian          = gradCons;
        funcs.jacobianstructure = @(x) sparse(cell2mat(consSparsity));
        %funcs.iterfunc         = @callback;
        
        % Set the IPOPT options.
        options.ipopt.hessian_approximation = 'limited-memory';
        options.ipopt.tol         = 1e-7;
        options.ipopt.max_iter    = 100;
        
        % Run IPOPT.
        [x info] = ipopt(x0,funcs,options);
        
        
    otherwise
        
        error('algorithm must be either remso, ipopt, or snopt')
        
        
end