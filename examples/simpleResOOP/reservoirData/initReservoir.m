function [reservoirP] = initReservoir( eclipseFile,varargin)
%
%Initialize a reservoir structure for REMSO.  A reservoir structure must contain:
%
%
%reservoirP.rock = rock;
%reservoirP.fluid = fluid;
%reservoirP.schedule = schedule;
%reservoirP.G = G;
%reservoirP.state = state;
%reservoirP.system = system;
%
%
% i.e.  the parameters to runScheduleADI.
%
%
%
require 

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


%% Compute constants
% Once we are happy with the grid and rock setup, we compute
% transmissibilities. For this we first need the centroids.
G = computeGeometry(G);

%% Set up reservoir
% We turn on gravity and set up reservoir and scaling factors.
gravity on

state = initResSol(G, deck.SOLUTION.EQUIL(2), [.15, .85]);

%{
system = initADISystem({'Oil', 'Water'}, G, rock, fluid);
%MRST-2014a eqsfiOW provide wrong jacobians WRT the wells
system.getEquations = @eqsfiOWExplicitWells;

system.well.allowControlSwitching = false;
system.well.allowCrossFlow = true;
system.well.allowWellSignChange = true;
system.well.cdpCalc = 'none';
%}


model = selectModelFromDeck(G, rock, fluid, deck);
schedule = convertDeckScheduleToMRST(G, model, rock, deck);
schedule.time = 0; % not set from the eclipse file


%{
mrstVerbose on
timer = tic;
[wellSols rSolOut] = simulateScheduleAD(state, model, schedule,'OutputMinisteps', true);
toc(timer);


plotRunScheduleSolution( state,rSolOut,schedule )

save forwardRun


%}



reservoirP.model = model;
reservoirP.schedule = schedule;
reservoirP.state = state;


end

