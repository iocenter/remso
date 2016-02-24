function obj = networkJointNPVConstraints(forwardStates,schedule,p, nCells, netSol, freqScale, pressureScale, flowScale, numStages, qlMin, qlMax, pScale, varargin)
% Compute net present value of a schedule with well solutions
% Inspired in NPVOW
% This function only changes the inputs, and have additional options

opt     = struct('OilPrice',             1.0 , ...
                'WaterProductionCost',  0.1 , ...
                'WaterInjectionCost',   0.1 , ...
                'DiscountFactor',       0.1 , ...
                'ComputePartials',      false, ...
                'tStep' ,               [],...
                'scale',                1    ,...
                'leftSeed',[],...
                'sign',1, ...        
                'turnoffPumps', false, ...
                'dpFunction', @simpleDp, ...
                'extremePoints', []);

opt     = merge_options(opt, varargin{:});

ro  = opt.OilPrice            / stb;
rw  = opt.WaterProductionCost / stb;
ri  = opt.WaterInjectionCost  / stb;
d   = opt.DiscountFactor;


wellSols = cellfun(@(x)x.wellSol,forwardStates,'UniformOutput',false);

% pressure and saturaton vectors just used for place-holding
pressure  = zeros(nCells, 1);
sW = zeros(nCells, 1);

dts   = schedule.step.val;

tSteps = opt.tStep;
if isempty(tSteps) %do all
    time = 0;
    numSteps = numel(dts);
    tSteps = (1:numSteps)';
else
    time = sum(dts(1:(opt.tStep-1)));
    numSteps = 1;
    dts = dts(opt.tStep);
end


%%%%%%%%%%%%%%%%%%%%%
%% Running network  %
%%%%%%%%%%%%%%%%%%%%%%
lastStep   = numel(forwardStates);     
wellSol = wellSols{lastStep};

netSol = runNetwork(netSol, wellSol, forwardStates{lastStep}, p, pScale, 'ComputePartials', opt.ComputePartials, 'turnoffPumps', opt.turnoffPumps, 'dpFunction', opt.dpFunction);   % running the network

vw = getVertex(netSol, netSol.VwProd);
ew = getEdge(netSol, vertcat(vw.Eout));
    
qf = vertcat(ew.qoE) + vertcat(ew.qwE);  % flows in the pumps   

%%%%%%%%%%%%%%%%%%%%
%% Equipment Dp   %%
%%%%%%%%%%%%%%%%%%%%
dpf = getChokesDp(netSol); % dp in the pumps

if  ~isempty(opt.extremePoints)
    qfWells = abs(qf)./(meter^3/day);
    dpPumps = abs(dpf)./barsa;
    
    qminFmin = cell2mat(opt.extremePoints(1));
    qminFmax = cell2mat(opt.extremePoints(2));
    qmaxFmin = cell2mat(opt.extremePoints(3));
    qmaxFmax = cell2mat(opt.extremePoints(4));
    
    [m1, n1] =  linearCoefficients(qminFmin, qmaxFmin);
    [m2, n2]=  linearCoefficients(qmaxFmin, qmaxFmax);
    [m3, n3] =  linearCoefficients(qminFmax, qmaxFmax);
    [m4, n4] =  linearCoefficients(qminFmin, qminFmax);
    
    c1 = dpPumps -(m1.*qfWells + n1);
    c2 = (m3.*qfWells + n3) - dpPumps;
    c3 = ((dpPumps - n2)./m2) - qfWells;    
    c4 = qfWells - ((dpPumps - n4)./m4);    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%    pump efficiency extreme points    %%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    rf = 0.75; % cost
    qmidFmin = (qminFmin(:,1) + qmaxFmin(:,1) )./2;
    qmidFmax = (qminFmax(:,1)  + qmaxFmax(:,1))./2;
    
    dpmFmin = m1.*qmidFmin + n1;
    dpmFmax = m3.*qmidFmax + n3;
    
    [mEfficiency, nEfficiency] = linearCoefficients([qmidFmin,dpmFmin], [qmidFmax,dpmFmax]);    
    penaltyEfficiency =  spones(ones(1, numel(freqScale)))*rf*(dpPumps -(mEfficiency.*qfWells + nEfficiency)).^2;     
    
    obj = cell(1,numSteps);
    for step = 1:numSteps
        sol = wellSols{tSteps(step)};
        qWs  = vertcat(sol.qWs);
        qOs  = vertcat(sol.qOs);
        injInx  = (vertcat(sol.sign) > 0);

        nW  = numel(qWs);
        pBHP = zeros(nW, 1); %place-holder

        if opt.ComputePartials
            [~, ~, qWs, qOs, ~,p] = initVariablesADI(pressure, sW, qWs, qOs, pBHP,p);
        end

        dt = dts(step);
        time = time + dt;

        prodInx = ~injInx;               
        
        objNPV = opt.scale*opt.sign*( dt*(1+d)^(-time/year) )*...
            spones(ones(1, nW))*( (-ro*prodInx).*qOs ...
            +(rw*prodInx - ri*injInx).*qWs ) + ...
             opt.scale*penaltyEfficiency;              

        if step<numSteps
            obj{step} = [c1.*0; c2.*0; c3*0; c4*0; dpf*0;  objNPV.*0];
        else
            obj{step} = [ [c3*(meter^3/day); c4*(meter^3/day)]./flowScale; c1*barsa./pressureScale; c2*barsa./pressureScale;  dpf./pressureScale; objNPV];
        end

        if opt.ComputePartials && ~(size(opt.leftSeed,2)==0)
            obj{step}.jac = cellfun(@(x)opt.leftSeed*x,obj{step}.jac,'UniformOutput',false);
        end
    end    
else 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Equipment Frequency   %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    rf = 0.1;
    fref = 60; %Hz 

    freq = 0./qf; % initialize frequency vector

    inletStr = vertcat(ew.stream);
    mixtureDen = (vertcat(inletStr.oil_dens) + vertcat(inletStr.water_dens))./2;  % density of the mixture

    dhf= pump_dh(dpf, mixtureDen); % dh in the pumps


    cond_dhf = dhf < 0;
    if any(cond_dhf)
        freq(cond_dhf) = pump_eq_system_explicit(qf(cond_dhf), dhf(cond_dhf), 60, numStages(cond_dhf));  % solves a system of equations to obtain frequency, flow and dh at 60Hz
    end

    if isa(freq, 'ADI')
        freq.val = real(freq.val);
        freq.jac = cellfun(@(w) real(w) ,freq.jac, 'UniformOutput', false);
    else
        freq = real(freq);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Flow Constraint in Equipment  %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    qpump_min = pump_rate(freq, qlMin, 60);
    pump_min = (-qf - qpump_min);


    qpump_max = pump_rate(freq, qlMax, 60);   
    pump_max = (qpump_max + qf);


    obj = cell(1,numSteps);
    for step = 1:numSteps
        sol = wellSols{tSteps(step)};
        qWs  = vertcat(sol.qWs);
        qOs  = vertcat(sol.qOs);
        injInx  = (vertcat(sol.sign) > 0);

        nW  = numel(qWs);
        pBHP = zeros(nW, 1); %place-holder

        if opt.ComputePartials
            [~, ~, qWs, qOs, ~,p] = initVariablesADI(pressure, sW, qWs, qOs, pBHP,p);
        end

        dt = dts(step);
        time = time + dt;

        prodInx = ~injInx;

%         objFreq = spones(ones(1, numel(freqScale)))*(rf*((freq-fref).^2));
% 
%         objNPV = opt.scale*opt.sign*( dt*(1+d)^(-time/year) )*...
%             spones(ones(1, nW))*( (-ro*prodInx).*qOs ...
%             +(rw*prodInx - ri*injInx).*qWs ) + ...
%             opt.scale*objFreq;  

        objNPV = opt.scale*opt.sign*( dt*(1+d)^(-time/year) )*...
            spones(ones(1, nW))*( (-ro*prodInx).*qOs ...
            +(rw*prodInx - ri*injInx).*qWs );
           

        if step<numSteps
            obj{step} = [pump_min.*0; pump_max.*0; freq*0./freqScale ; dpf*0;  objNPV];
        else
            obj{step} = [ [pump_min; pump_max]./flowScale; freq./freqScale; dpf./pressureScale; objNPV];
        end

        if opt.ComputePartials && ~(size(opt.leftSeed,2)==0)
            obj{step}.jac = cellfun(@(x)opt.leftSeed*x,obj{step}.jac,'UniformOutput',false);
        end
    end
    
end

end



