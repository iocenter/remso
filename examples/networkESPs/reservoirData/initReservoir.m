function [reservoirP,units] = initReservoir( eclipseFile,varargin)
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
if(nargin < 1)
   eclipseFile = 'SIMPLE10x5x10.txt';
end

% input/output-files ------------------------------------------------------
cntrNm  = 'SIMPLE10x5x10_CONTROLS.TXT';
outNm   = 'SIMPLE10x5x10_RES.TXT';
gradNm  = 'SIMPLE10x5x10_GRAD.TXT';
% -------------------------------------------------------------------------


opt = struct('Verbose', false, ...
             'netWells', []);
         
opt = merge_options(opt, varargin{:});
verbose = opt.Verbose;


current_dir = fileparts(mfilename('fullpath'));
fn    = fullfile(current_dir, eclipseFile);
deck = readEclipseDeck(fn);


if isfield(deck.RUNSPEC, 'METRIC')
    units = 'METRIC';
elseif isfield(deck.RUNSPEC,  'FIELD')
    units = 'FIELD';
elseif isfield(deck.RUNSPEC, 'LAB')
    units = 'LAB';
else 
    error('Please specify explicitly the units in the Eclipse input file')
end


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


system = initADISystem({'Oil', 'Water'}, G, rock, fluid, 'netWells', opt.netWells);
%MRST-2014a eqsfiOW provide wrong jacobians WRT the wells
%system.getEquations = @eqsfiOWExplicitWells;

system.well.allowControlSwitching = false;
system.well.allowCrossFlow = true;
system.well.allowWellSignChange = true;
system.well.cdpCalc = 'none';

%if control input is given, edit schedule:
if ~isempty(dir(cntrNm));
    fid = fopen(cntrNm);
    u   = fscanf(fid, '%g');
    numWells      = numel(u)/numContrSteps;
    vals     = mat2cell(u(:), numWells*ones(numContrSteps, 1), 1);
    schedule = updateSchedule(schedule, vals);
    fclose(fid);
end


[schedule] = eclipseSchedule2mrstSchedule(schedule,G,rock);
[ schedule ] = relaxLimsInSchedule( schedule);

% % % % 
%schedule.step.val =  schedule.step.val(1:10);
%schedule.step.control = schedule.step.control(1:10);
%schedule.control = schedule.control(1:schedule.step.control(10));

%{
mrstVerbose on
timer = tic;
[wellSols rSolOut] = runScheduleADI(state, G, rock, system, schedule);
toc(timer);



[qWs, qOs, qGs, bhp] = wellSolToVector(wellSols);


save forwardRun


%}

reservoirP.rock = rock;
reservoirP.fluid = fluid;
reservoirP.schedule = schedule;
reservoirP.G = G;
reservoirP.state = state;
%reservoirP.scalFacs = scalFacs;
reservoirP.system = system;

end
