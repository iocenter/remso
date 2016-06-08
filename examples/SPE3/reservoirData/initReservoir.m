function [reservoirP,units] = initReservoir()


%% SPE3 case using fully implicit black oil solver
% This <http://dx.doi.org/10.2118/12278-PA third SPE comparative solution project>
% consists of a gas injection problem in a small (9x9x4) reservoir. The problem is
% originally set up to be solved using a compositional solver. Using PVTi, the physical
% datas have been processed to obtain an equivalent blackoil problem and the resulting
% blackoil parameters are provided in the file "SPE3.DATA". The data set we provide is a
% modified version of input files belonging to the
% <http://www.ntnu.edu/studies/courses/TPG4535 course in reservoir engineering and
% petrophysics> at NTNU (Trondheim, Norway) and available at
% <http://www.ipt.ntnu.no/~kleppe/pub/SPE-COMPARATIVE/ECLIPSE_DATA/>. The oil can vaporize
% but the gas cannot dissolve in oil so that the gas/oil ratio remains equal to zero
% during the whole simulation. The data are


current_dir = fileparts(mfilename('fullpath'));
fn    = fullfile(current_dir,'SPE3.DATA');

deck = readEclipseDeck(fn);

if isfield(deck.RUNSPEC, 'FIELD');
    units = 'FIELD';
elseif isfield(deck.RUNSPEC, 'LAB');
    units = 'LAB';    
elseif isfield(deck.RUNSPEC, 'SI');
	units = 'SI';
else
    units = 'METRIC'; % eclipse default
end
% The deck is given in field units, MRST uses consistently the metric system.
deck = convertDeckUnits(deck);

G = initEclipseGrid(deck);
if isfield('ACTNUM',deck.GRID)
    G     = removeCells(G, ~deck.GRID.ACTNUM);
end
G = computeGeometry(G);

rock  = initEclipseRock(deck);
rock  = compressRock(rock, G.cells.indexMap);

% Create a special ADI fluid which can produce differentiated fluid properties.
fluid = initDeckADIFluid(deck);

% gravity on
assert((norm(gravity) - 9.80665) == 0);  %% check that the gravity is on by default


%% Set up initial state
% The initial state corresponds to an equilibrium state between gravitational and
% capillary forces. It is now loaded from deck.

p0  = deck.SOLUTION.PRESSURE;
sw0 = deck.SOLUTION.SWAT;
sg0 = deck.SOLUTION.SGAS;
s0  = [sw0, 1-sw0-sg0, sg0];
rv0 = deck.SOLUTION.RV;
rs0 = 0;


if ~isfield('rsSat',fluid)
    fluid.rsSat = rs0;
end

state = struct('s', s0, 'rs', rs0, 'rv', rv0, 'pressure', p0);   


clear k p0 s0 rv0 rs0




schedule = deck.SCHEDULE;

[schedule] = eclipseSchedule2mrstSchedule(schedule,G,rock);
schedule  = relaxLimsInSchedule( schedule);
% The schedule given in deck provides negative flow without the limits!
schedule.control(2).W(1).val = schedule.control(1).W(1).val;
if ~isfield(schedule,'time')
    schedule.time = 0;
end

times = [1 20 70.25 91.25 91.25 91.25]'*day;  %% a year

scheduleF = schedule;
scheduleF.step.val = times;
scheduleF.step.control = ones(numel(times),1);
scheduleF.control = scheduleF.control(1);

scheduleL = schedule;
scheduleL.step.val = times;
scheduleL.step.control = ones(numel(times),1);
scheduleL.control = scheduleL.control(2);
schedules = [repmat(scheduleF,10,1);repmat(scheduleL,10,1)];

schedule = mergeSchedules(schedules);


system.well.allowControlSwitching = false;
system.well.allowWellSignChange = true;
system.well.allowCrossFlow = true;




system = initADISystem(deck, G, rock, fluid, 'cpr', false);
system.nonlinear.itLinearSolver = false;

% correct state!
[state] = updateStateVO([], state, repmat({zeros(size(state.pressure))},3,1), fluid, system,'updateWellSol',false);

system.well.cdpCalc = 'none';

reservoirP.state = state;
reservoirP.G = G;
reservoirP.rock = rock;
reservoirP.system = system;
reservoirP.schedule = schedule;
reservoirP.fluid = fluid;


%{
mrstVerbose on
timer = tic;
[wellSols, states, iter] = runScheduleADI(state, G, rock, system, schedule);
toc(timer) 


[cellfun(@(si)max(si.pressure)/barsa,states),cellfun(@(si)min(si.pressure)/barsa,states)]
[cellfun(@(si)max(si.s(:,1)),states),cellfun(@(si)min(si.s(:,1)),states)]

[qWs, qOs, qGs, bhp] = wellSolToVector(wellSols);

[max(qWs*day);min(qWs*day)]
[max(qOs*day);min(qOs*day)]
[max(qGs*day);min(qGs*day)]
[max(bhp/barsa);min(bhp/barsa)]


w = functions(fluid.rvSat);
pmax = w.workspace{1}.pvtg{1}.key(end);
pmin = w.workspace{1}.pvtg{1}.key(1);

px = (0:100)*(pmax-pmin)/100 + pmin;

plot(px,fluid.rvSat(px));

rvSatx = [ (0:99)*fluid.rvSat(pmin)/100 , fluid.rvSat(px) ]

rh = @(rv) fluid.rhoGS./(fluid.rhoGS+fluid.rhoOS*rv);

plot(rvSatx,rh(rvSatx))

mrstVerbose off
%}

end
