function obj = networkJointNPVConstraints(forwardStates,schedule, nCells, netSol, freqScale, pressureScale, flowScale, numStages, fref, qlMin, qlMax, varargin)
%%
% networkJointNPVConstraints computes net present value and the network algebraic variables
% Inspired in NPVOW
% This function only changes the inputs, and have additional options
%
% SYNOPSIS:
%  [u,x,v,f,xd,M,simVars] = networkJointNPVConstraints(forwardStates,schedule, nCells, netSol, freqScale, pressureScale, flowScale, numStages, fref, qlMin, qlMax, ...)
% PARAMETERS:
%   forwardStates - Reservoir states at a given prediction interval
%   schedule - Schedule for the given prediction interval
%   freqScale, pressureScale, flowScale - Scale objects for pump frequency, pressure and flow rates 
%   numStages,- number of stages of the ESP
%   fref - ESP reference frequency
%   qlMin, qlMax - ESP flow boundaries for the reference frequency
% RETURNS:
%
%  obj - vector with values and gradients of objective and network variables
%
%{

Copyright 2015-2018, Thiago Lima Silva

REMSO is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

REMSO is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with REMSO.  If not, see <http://www.gnu.org/licenses/>.

%}
opt     = struct('OilPrice',             1.0 , ...
                'WaterProductionCost',  0.1 , ...
                'WaterInjectionCost',   0.1 , ...
                'DiscountFactor',       0.1 , ...
                'ComputePartials',      false, ...
                'tStep' ,               [],...
                'scale',                1    ,...
                'leftSeed',[],...
                'sign',1, ...                        
                'dpFunction', @dpBeggsBrillJDJ, ...
                'forwardGradient',true,...
                'finiteDiff', true,...                
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

netSol = runNetwork(netSol, wellSol, forwardStates{lastStep}, 'ComputePartials', opt.ComputePartials, 'dpFunction', opt.dpFunction, 'forwardGradient', opt.forwardGradient,'finiteDiff', opt.finiteDiff);   % running the network

if  ~isempty(opt.extremePoints) %% linear approximation of pump constraints
    qf = netSol.qo(netSol.Eeqp) + netSol.qw(netSol.Eeqp);  % flows in the pumps    
    
    %%%%%%%%%%%%%%%%%%%%
    %% Equipment Dp   %%
    %%%%%%%%%%%%%%%%%%%%
    equip = getEdge(netSol,netSol.Eeqp);
    vin = vertcat(equip.vin);
    vout = vertcat(equip.vout);
    
    dpf = netSol.pV(vin)-netSol.pV(vout); % dp in the pumps    
    
    qfWells = abs(qf)./(meter^3/day); %% TODO: check direction of the flow
    dpPumps = abs(dpf)./barsa;        %% TODO: check is pressure drop is positive or negative for the pump
    
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
%     rf = 0.75; % cost
    rf = 0; % no cost function
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
            [~, ~, qWs, qOs, ~] = initVariablesADI(pressure, sW, qWs, qOs, pBHP);
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
else  %% original nonlinear pump constraints    
    [freq, qf, qpump_min, qpump_max, ~, ~] = nonlinearPumpConstraints(netSol, fref, numStages, qlMin, qlMax);
    
    pump_min = (-qf - qpump_min);
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
            [~, ~, qWs, qOs, ~] = initVariablesADI(pressure, sW, qWs, qOs, pBHP);
        end

        dt = dts(step);
        time = time + dt;

        prodInx = ~injInx;

        objNPV = opt.scale*opt.sign*( dt*(1+d)^(-time/year) )*...
            spones(ones(1, nW))*( (-ro*prodInx).*qOs ...
            +(rw*prodInx - ri*injInx).*qWs );
            
        

        if step<numSteps
            obj{step} = [pump_min.*0; pump_max.*0; freq*0./freqScale ; objNPV];
        else
            obj{step} = [ [pump_min; pump_max]./flowScale; freq./freqScale; objNPV];
        end

        if opt.ComputePartials && ~(size(opt.leftSeed,2)==0)
            obj{step}.jac = cellfun(@(x)opt.leftSeed*x,obj{step}.jac,'UniformOutput',false);
        end
    end
    
end

end



