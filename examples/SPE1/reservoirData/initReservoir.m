function [reservoirP] = initReservoir( eclipseFile,varargin)


if nargin < 1
    eclipseFile = 'odeh_adi.data';
end

opt = struct('Verbose', false);
         
opt = merge_options(opt, varargin{:});
verbose = opt.Verbose;

%% SPE1 case 
% This <http://dx.doi.org/10.2118/9723-PA first comparative solution project> consists of
% a gas injection problem in a small ($10\times10\times3$) reservoir with a single
% producer and a single injector. It is set up to be solved using a black-oil model. The
% data set we provide is a modified version of input files belonging to the
% <http://www.ntnu.edu/studies/courses/TPG4535 course in reservoir engineering and
% petrophysics> at NTNU (Trondheim, Norway) and available at
% <http://www.ipt.ntnu.no/~kleppe/pub/SPE-COMPARATIVE/ECLIPSE_DATA/>. The results are
% compared to the output from a major commercial reservoir simulator (Eclipse 100).


%% Read input files
% The input files follow Eclipse format. MRST contains a dedicated module which can handle
% standard Eclipse keywords.

mrstVerbose off


deck = readEclipseDeck(eclipseFile);

% The deck is given in field units, MRST uses metric.
deck = convertDeckUnits(deck,'verbose',verbose);

G = initEclipseGrid(deck);
G = computeGeometry(G);

rock  = initEclipseRock(deck);
rock  = compressRock(rock, G.cells.indexMap);

% Create a special ADI fluid which can produce differentiated fluid properties.
fluid = initDeckADIFluid(deck);

% The case includes gravity
gravity on

%% Setup initial state
% The initial state is a pressure field that is constant in each layer, a
% uniform mixture of water (Sw=0.12) and oil (So=0.88) with no initial free
% gas (Sg=0.0) and a constant dissolved gas/oil ratio (|Rs|) throughout the
% model. The pressure and Rs values are derived through external means.

[k, k, k] = gridLogicalIndices(G);

p0    = [ 329.7832774859256 ; ...  % Top layer
          330.2313357125603 ; ...  % Middle layer
          330.9483500720813 ];     % Bottom layer

p0    = convertFrom(p0(k), barsa);
s0    = repmat([ 0.12, 0.88, 0.0 ], [G.cells.num, 1]);
rs0   = repmat( 226.1966570852417 , [G.cells.num, 1]);
rv0   = 0; % dry gas

if ~isfield('rvSat',fluid)
    fluid.rvSat = rv0;
end


state = struct('s', s0, 'rs', rs0, 'rv', rv0, 'pressure', p0);
clear k p0 s0 rs0;



%% Initialize schedule and system before solving for all timesteps
% We extract the schedule from the read deck and create a ADI system for our problem. The
% system autodetects a black oil problem and sets up default values for the various
% options. The only thing we change is that we disable the CPR preconditioner as the
% problem is too small to benefit from preconditioning: The overhead required for the
% preconditioner is bigger than the benefits offered by a specialized solver.
%
% During some time steps (67 and 91) the Newton iterations oscillate. The
% solver detects this, and dampens or relaxes the step length when this
% behavior is observed.
%
% To see detailed convergence analysis during each time step, set verbose
% to on by using: |mrstVerbose on|

schedule = deck.SCHEDULE;

if ~isfield(schedule,'time')
    schedule.time = 0;
end

system = initADISystem(deck, G, rock, fluid, 'cpr', false);

system.well.allowControlSwitching = false;
system.well.allowCrossFlow = true;
system.well.allowWellSignChange = false;
system.well.cdpCalc = 'none';
system.nonlinear.linesearch = true;

system.stepOptions.dpMax = 25*barsa;
system.stepOptions.drsMax = 25;
system.stepOptions.dsMax = 0.1;
system.stepOptions.minp = 1000*psia;

[schedule] = eclipseSchedule2mrstSchedule(schedule,G,rock);


[ schedule ] = relaxLimsInSchedule( schedule);

%schedule.control.W(1).val =  100*(1000 * ft^3/day);    
%schedule.control.W(2).val = -1000*stb/day;

stepSize = 5;
schedule.step.val = stepSize*day*ones(ceil(1200/stepSize),1);
schedule.step.control = ones(ceil(1200/stepSize),1);



% cut the existing schedule in 10
nIntervals = 10;
intervals = ceil((1/nIntervals:1/nIntervals:1)*numel(schedule.step.val));

schedules = multipleSchedules(schedule, intervals );

for k = 1:numel(schedules)
    schedules(k).step.val = [schedules(k).step.val];
    schedules(k).step.control = [schedules(k).step.control];
end
schedule = mergeSchedules(schedules);


reservoirP.rock = rock;
reservoirP.fluid = fluid;
reservoirP.schedule = schedule;
reservoirP.G = G;
reservoirP.state = state;
reservoirP.system = system;

end

