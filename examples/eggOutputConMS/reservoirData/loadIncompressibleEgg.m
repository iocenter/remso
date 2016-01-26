function [ reservoirP ] = loadIncompressibleEgg( eggDir )


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
% Gijs van Essen	SPE 124332 Model properties (Report parameters)
fluid     = initCoreyFluid('mu' ,  [   1,   5]*centi*poise     , ...
                            'rho',  [1000, 900]*kilogram/meter^3, ...
                            'n'  ,  [   3,   4]                 , ...
                            'sr' ,  [ 0.2,   0.1]                 , ...
                            'kwm',  [   0.749,0.8]);
% observe that pc (capillary pressure) is not supported when computing gradients
                        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%
%%%     Reading rock properties (permeability and porosity fields)
%%%
      
rock  = initEclipseRock(deck);
rock  = compressRock(rock, G.cells.indexMap);

gravity off 
assert((norm(gravity) - 0) == 0);  %% gravity will be diregarded an
   
% Get schedule
fprintf('Overwriting schedule given on deck!\n');
schedule = deck.SCHEDULE;

schedule = eclipseSchedule2mrstSchedule(schedule,G,rock);
W = schedule.control(1).W;
W = assembleWellSystem(G, W, 'Type', 'comp_hybrid');

nYears = 10;
ystep = [1,4,10,15,30*ones(1,11)]'*day;
simulationControlSteps = rldecode((1:nYears)', numel(ystep)*ones(nYears,1));
simulationTimeSteps  = repmat(ystep, nYears, 1);
simulationTimeSteps = cumsum(simulationTimeSteps);

%% Initialize and construct the linear system
S            = computeMimeticIP(G, rock,'Type', 'tpfa','InnerProduct', 'ip_tpf','Verbose', true); 


% Change of pressure initializaton!!
po(1) = 400*barsa;

rSol         = initResSol(G, po(1),0.1);
rSol.wellSol = initWellSol(W, po(1)); 

if strcmp(S.type, 'hybrid')
   solver = 'hybrid';
else
   solver = 'mixed';
end
rSol = solveIncompFlowLocal(rSol, G, S, fluid, 'wells',W,'Solver',solver);

schedule = initSchedule(W, 'TimeSteps', simulationTimeSteps, 'Verbose', false);


%{
[simRes,reports] = runSchedule(rSol, G, S, W, rock, fluid,schedule,'Verbose',false ,'VerboseLevel',2,'Verbose',true);


 [qw,qo,bhp,t0,tf] = perfVariables(W, fluid, simRes);

t = [t0,tf];
t = reshape(t',numel(t),1);
qwP = cell2mat(cellfun(@(q)[q,q]',qw,'UniformOutput',false));
qoP = cell2mat(cellfun(@(q)[q,q]',qo,'UniformOutput',false));
bhpP = cell2mat(cellfun(@(q)[q,q]',bhp,'UniformOutput',false));

figure(1);plot(t/day,qwP*day)
figure(2);plot(t/day,qoP*day)
figure(3);plot(t/day,bhpP/barsa)


%}






reservoirP.rock = rock;
reservoirP.fluid = fluid;
reservoirP.schedule = schedule;
reservoirP.G = G;
reservoirP.state = rSol;
reservoirP.S = S;
reservoirP.W = W;
reservoirP.simulationControlSteps = simulationControlSteps;







end

