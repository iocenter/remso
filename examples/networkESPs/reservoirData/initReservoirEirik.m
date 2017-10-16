function [reservoirP,units] = initReservoir(varargin)
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
% input/output-files ------------------------------------------------------
caseNm  = 'SIMPLE10x5x10.txt';
cntrNm  = 'SIMPLE10x5x10_CONTROLS.TXT';
outNm   = 'SIMPLE10x5x10_RES.TXT';
gradNm  = 'SIMPLE10x5x10_GRAD.TXT';
% -------------------------------------------------------------------------

% Enable this to get convergence reports when solving schedules
mrstVerbose false

% check if initialized model already exists
[pth, nm, ext] = fileparts(caseNm);
% modelNm = fullfile(pth,[nm,'.mat']);
% doInitialize = true;
% if doInitialize
    deck = readEclipseDeck(caseNm);
    % Convert to MRST units (SI)
    deck = convertDeckUnits(deck);
    % Create grid
    G = initEclipseGrid(deck);
    % Set up the rock structure
    rock  = initEclipseRock(deck);
    rock  = compressRock(rock, G.cells.indexMap);
    % Create fluid
    fluid = initDeckADIFluid(deck);
    % Get schedule
    schedule = deck.SCHEDULE;

    %% Compute constants
    G = computeGeometry(G);

    %% Set up reservoir
    % We turn on gravity and set up reservoir and scaling factors.
    gravity on
    state = initStateADI(G, fluid, deck);
    system = initADISystem({'Oil', 'Water'}, G, rock, fluid);

    system.well.cdpCalc = 'none';
    
    % save initialized model
%     save(modelNm, 'G', 'state', 'rock', 'fluid', 'schedule', 'system');
% else
%     load(modelNm, 'G', 'state', 'rock', 'fluid', 'schedule', 'system');
% end

%if control input is given, edit schedule:
numContrSteps = numel(schedule.control);
if ~isempty(dir(cntrNm));
    fid = fopen(cntrNm);
    u   = fscanf(fid, '%g');
    numWells      = numel(u)/numContrSteps;
    vals     = mat2cell(u(:), numWells*ones(numContrSteps, 1), 1);
    schedule = updateSchedule(schedule, vals);
    fclose(fid);
end

[schedule] = eclipseSchedule2mrstSchedule(schedule,G,rock);
schedule.time = 0;


schedule.step.val =  schedule.step.val(1:10);
schedule.step.control = schedule.step.control(1:10);
schedule.control = schedule.control(1:schedule.step.control(10));

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

