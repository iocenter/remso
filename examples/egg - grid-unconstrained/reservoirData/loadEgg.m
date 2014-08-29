function [ reservoirP ] = loadEgg( eggDir )


fn    = fullfile(eggDir, ...
    'Egg_Model_ECL.DATA');

if ~exist(fn, 'file'),
   error('Egg model data is not available.')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% 
%%%     Reading the input deck
%%%
deck  = readEclipseDeck(fn);      
deck = convertDeckUnits(deck);  % Convert to MRST units (SI)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%
%%%     Reading grid structure
%%%
G     = initEclipseGrid(deck); 
G     = removeCells(G, ~deck.GRID.ACTNUM);
G     = computeGeometry(G);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%
%%%     Defining fluid properties
%%%
fluid = initDeckADIFluid(deck);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%
%%%     Reading rock properties (permeability and porosity fields)
%%%
      
rock  = initEclipseRock(deck);
rock  = compressRock(rock, G.cells.indexMap);

% gravity on 
assert((norm(gravity) - 9.80665) == 0);  %% check that the gravity is on by default
   
% Get schedule
fprintf('Overwriting schedule given on deck!\n');
schedule = deck.SCHEDULE;



% make 10-year schedule with 10-control steps
ystep = [1,4,10,15,30*ones(1,11)]'*day;
s.step.val     = repmat(ystep, 10, 1);
s.step.control = rldecode((1:10)', numel(ystep)*ones(10,1));


schedule.control = repmat(schedule.control,10,1);
schedule.step = s.step;
schedule.time = 0;


% Introduce wells
W = processWells(G, rock, deck.SCHEDULE.control(1)); % Introduce wells


% Change of pressure initializaton!!
po(1) = 400*barsa;


rSol         = initResSol(G, po(1),0.1);

%%fix bug
if size(rSol.s,2) == 1
    sat = zeros(size(rSol.s,1),2);
    sat(:,1) = rSol.s;
    sat(:,2) = 1 - sat(:,1);
    rSol.s = sat;
end

 
system = initADISystem({'Water', 'Oil'}, G, rock, fluid);


% system setings:
system.nonlinear.cpr = true;
% use new cpr based on dnr
system.nonlinear.cprBlockInvert = false;
% convergence is overall better for quite strict limits on update 
system.stepOptions.drsMax = .5;
system.stepOptions.dpMax  = .5;
system.stepOptions.dsMax  = .2;
% gmres tol needs to be quite strict
system.nonlinear.cprRelTol = 1e-3;
system.pscale = 1/(400*barsa);   
% use direct solver instead !!!
system.nonlinear.itLinearSolver = false;


system.well.allowControlSwitching = false;
system.well.allowCrossFlow = true;
system.well.allowWellSignChange = true;

[schedule] = eclipseSchedule2mrstSchedule(schedule,G,rock);


[ schedule ] = relaxLimsInSchedule( schedule);



%{
mrstVerbose on
timer = tic;
[wellSols rSolOut] = runScheduleADI(rSol, G, rock, system, schedule);
toc(timer);
%}



reservoirP.rock = rock;
reservoirP.fluid = fluid;
reservoirP.schedule = schedule;
reservoirP.G = G;
reservoirP.state = rSol;
reservoirP.system = system;









end

