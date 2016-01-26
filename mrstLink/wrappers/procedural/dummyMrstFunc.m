function obj = dummyMrstFunc(state,wellSol,schedule,varargin)
%
% This function shows the most general structure of a mrst point function
% considered in the Remso algorithm.  It also shows the structure of the
% Jacobian.
%
% state refers to the states at the end of the step.
% wellsol is the wellSol during the simulation period
% schedule gives the control values
%

opt     = struct('ComputePartials', false);
opt     = merge_options(opt, varargin{:});




% pressure and saturaton vectors just used for place-holding
p  = state.pressure;
sW = state.s(:,1);

pBHP = vertcat(wellSol.bhp);
schVal = schedule2Controls(schedule); 
qWs  = vertcat(wellSol.qWs);
qOs  = vertcat(wellSol.qOs);

if opt.ComputePartials
    [p, sW, qWs, qOs, pBHP, schVal] = initVariablesADI(p, sW, qWs, qOs, pBHP,schVal);
end

% f can have several lines (outputs), for fun here we make the jacobian to
% be [1,2,3..]
nV = 0;
obj = 0;

n = numel(double(p));
obj = obj + sum((nV+1:nV+n)'.*p);
nV = nV + n;

n = numel(double(sW));
obj = obj + sum((nV+1:nV+n)'.*sW);
nV = nV + n;

n = numel(double(qWs));
obj = obj + sum((nV+1:nV+n)'.*qWs);
nV = nV + n;

n = numel(double(qOs));
obj = obj + sum((nV+1:nV+n)'.*qOs);
nV = nV + n;

n = numel(double(pBHP));
obj = obj + sum((nV+1:nV+n)'.*pBHP);
nV = nV + n;

n = numel(double(schVal));
obj = obj + sum((nV+1:nV+n)'.*schVal);
