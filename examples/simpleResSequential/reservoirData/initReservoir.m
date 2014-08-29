function [reservoirP] = initReservoir( eclipseFile,varargin)
%
%Initialize a reservoir structure for REMSO.  A reservoir structure must contain:
%
%
%reservoirP.rock = rock;
%reservoirP.fluid = fluid;
%reservoirP.schedule = schedule;
%reservoirP.G = G;
%reservoirP.T = T;
%reservoirP.state = state;
%reservoirP.system = system;
%
%
% i.e.  the parameters to runScheduleADI.
%
%
%
if(nargin < 1)
   eclipseFile = 'simple10x1x10.data';
end

opt = struct('Verbose', false);
         
opt = merge_options(opt, varargin{:});
verbose = opt.Verbose;


current_dir = fileparts(mfilename('fullpath'));
fn    = fullfile(current_dir, eclipseFile);
deck = readEclipseDeck(fn);

% Convert to MRST units (SI)
deck = convertDeckUnits(deck,'verbose',verbose);

% Create grid
G = initEclipseGrid(deck);

% Set up the rock structure
rock  = initEclipseRock(deck);
rock  = compressRock(rock, G.cells.indexMap);

% Create fluid
fluid = initDeckADIFluid(deck);

% Get schedule
schedule = deck.SCHEDULE;
schedule.time = 0; % not set from the eclipse file


%% Compute constants
% Once we are happy with the grid and rock setup, we compute
% transmissibilities. For this we first need the centroids.
G = computeGeometry(G);

%% Set up reservoir
% We turn on gravity and set up reservoir and scaling factors.
gravity on

state = initResSol(G, deck.SOLUTION.EQUIL(2), [.15, .85]);


system = initADISystem({'Oil', 'Water'}, G, rock, fluid);


system.well.allowControlSwitching = false;
system.well.allowCrossFlow = true;
system.well.allowWellSignChange = true;

[schedule] = eclipseSchedule2mrstSchedule(schedule,G,rock);


[ schedule ] = relaxLimsInSchedule( schedule);


reservoirP.rock = rock;
reservoirP.fluid = fluid;
reservoirP.schedule = schedule;
reservoirP.G = G;
reservoirP.state = state;
%reservoirP.scalFacs = scalFacs;
reservoirP.system = system;

end

