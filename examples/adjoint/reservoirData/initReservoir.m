function [ reservoirP ] = initReservoir( )
%INITRESERVOIR Summary of this function goes here
%   Detailed explanation goes here



% whether or not to show output
verbose = false;

% Define model ------------------------------------------------------------
nx = 21; ny = 21; nz = 1;
G = cartGrid([nx ny nz], [5*nx 5*ny 1*nz]);
G = computeGeometry(G);

c = G.cells.centroids;
rock.perm  = max(10*sin(c(:,2)/25+.5*cos(c(:,1)/25))-9, .01)*1000*milli*darcy;
rock.poro  = repmat(0.3, [G.cells.num, 1]);

fluid  = initCoreyFluid('mu' , [1, 5] .* centi*poise, ...
    'rho', [1014, 859].*kilogram/meter^3, ...
    'n'  , [2, 2], 'sr', [0, 0], 'kwm', [1, 1]);
fluid  = adjointFluidFields(fluid);


% Wells and initial rates -------------------------------------------------
radius = .1;
totTime = 500*day;
W = [];
% Injectors along left side:
nInj = 3; % > 1
pos  = (1 : (ny-1)/(nInj-1) : ny)';
posInj  = round(pos);
for k = 1:nInj
    nm = ['inj', num2str(k)];
    W = addWell(W, G, rock, 1+(posInj(k)-1)*nx, 'Type', 'bhp' , 'Val', 500*barsa, ...
        'Radius', radius, 'Name', nm, 'Comp_i', [1, 0], 'Sign', 1, 'InnerProduct', 'ip_tpf');
end
% Producers along right side:
nProd = 5; % >1
pos  = (1 : (ny-1)/(nProd-1) : ny)';
posProd  = round(pos);
for k = 1:nProd
    nm = ['prod', num2str(k)];
    W = addWell(W, G, rock, nx+(posProd(k)-1)*nx, 'Type', 'bhp' , 'Val', 150*barsa, ...
        'Radius', radius, 'Name', nm, 'Comp_i', [1, 0], 'Sign', -1, 'InnerProduct', 'ip_tpf');
end

% System components -------------------------------------------------------
S = computeMimeticIP(G, rock, 'Type', 'comp_hybrid', 'Verbose', verbose, ...
    'InnerProduct', 'ip_tpf');
W = assembleWellSystem(G, W, 'Type', 'comp_hybrid');

% Initialize --------------------------------------------------------------
state = initResSol(G, 0.0,0.02);
state.wellSol = initWellSol(W, 0);

% Objective function ------------------------------------------------------

% Initialize schedule and controls ----------------------------------------
numSteps = 10;

schedule = initSchedule(W, 'NumSteps', numSteps, 'TotalTime', ...
    totTime, 'Verbose', verbose);


reservoirP.state = state;
reservoirP.G = G;
reservoirP.S = S;
reservoirP.W = W;
reservoirP.rock = rock;
reservoirP.fluid = fluid;
reservoirP.schedule = schedule;

end

