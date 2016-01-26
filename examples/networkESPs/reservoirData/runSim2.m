function [] = runSim2

mrstModule add ad-fi
mrstModule add deckformat

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
modelNm = fullfile(pth,[nm,'.mat']);
doInitialize = isempty(dir(modelNm));
if doInitialize
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

    % save initialized model
    save(modelNm, 'G', 'state', 'rock', 'fluid', 'schedule', 'system');
else
    load(modelNm, 'G', 'state', 'rock', 'fluid', 'schedule', 'system');
end

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

% [schedule] = eclipseSchedule2mrstSchedule(schedule,G,rock);

%run simulation:
[wellSols, states] = runScheduleADI(state, G, rock, system, schedule);

%produce output
numWells = numel(wellSols{1});
numSteps = numel(wellSols);
wrats    = zeros(numWells, numSteps);
orats    = zeros(numWells, numSteps);
bhps     = zeros(numWells, numSteps);
for wn = 1:numWells
    wrats(wn,:) = cellfun(@(x)x(wn).qWs, wellSols)';
    orats(wn,:) = cellfun(@(x)x(wn).qOs, wellSols)';
    bhps(wn,:)  = cellfun(@(x)x(wn).pressure, wellSols)';
end

% take average over controlsteps
M = sparse((1:numSteps)', schedule.step.control, schedule.step.val, numSteps, numContrSteps);
wrats = (wrats*M)./(ones(numWells, 1)*sum(M));
orats = (orats*M)./(ones(numWells, 1)*sum(M));
bhps  = (bhps*M)./(ones(numWells, 1)*sum(M));

% write output
fid = fopen(outNm, 'w');
fprintf(fid, '%+12.6e %+12.6e %+12.6e\n',  [wrats(:) orats(:) bhps(:)]');
fclose(fid);

end

