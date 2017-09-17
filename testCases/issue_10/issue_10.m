function [ok] = issue_10()

    % REservoir Multiple Shooting Optimization.
    % REduced Multiple Shooting Optimization.

    % Make sure the workspace is clean before we start
    clc
    clear
    clear global

    % Required MRST modules
    mrstModule clear;
    mrstModule add deckformat ad-core ad-blackoil ad-props;

    here = fileparts(mfilename('fullpath'));
    if isempty(here)
        here = pwd();
    end
    warning('off','MATLAB:rmpath:DirNotFound')
    remso_root = fullfile(here,filesep,'..',filesep,'..',filesep);
    rmpath(genpath(remso_root));
    warning('on','MATLAB:rmpath:DirNotFound')

    
    % Include REMSO functionalities
    addpath(genpath(fullfile(remso_root,'mrstDerivated')));
    addpath(genpath(fullfile(remso_root,'mrstLink',filesep,'utils')));
    addpath(genpath(fullfile(remso_root,'mrstLink',filesep,'wrappers',filesep,'procedural')));
    addpath(genpath(fullfile(remso_root,'mrstLink',filesep,'wrappers',filesep,'OOP')));
    addpath(genpath(fullfile(remso_root,'optimization',filesep,'utils')));
    addpath(genpath(fullfile(remso_root,'optimization',filesep,'multipleShooting')));
    addpath(genpath(fullfile(remso_root,'optimization',filesep,'testFunctions')));
    addpath(genpath(fullfile(here,filesep,'reservoirData')));

    %% Initialize reservoir -  the Simple reservoir
    [reservoirP] = initReservoir( 'reallySimpleRes.data','Verbose',true);


    reservoirP.model.toleranceMB = 1e-12;
    reservoirP.model.toleranceCNV = 1e-7;
    reservoirP.model.toleranceWellBHP = barsa/1e-5;
    reservoirP.model.toleranceWellRate = 1/day/1e-5;  
    
    model = reservoirP.model;

    % do not display reservoir simulation information!
    mrstVerbose off;

    %% Multiple shooting problem set up
    totalPredictionSteps = numel(reservoirP.schedule.step.val);  % MS intervals

    % Schedule partition for each control period and for each simulated step
    lastControlSteps = findControlFinalSteps( reservoirP.schedule.step.control );
    controlSchedules = multipleSchedules(reservoirP.schedule,lastControlSteps);

    stepSchedules = multipleSchedules(reservoirP.schedule,1:totalPredictionSteps);


    %% Variables Scaling
    cellControlScales = schedules2CellControls(schedulesScaling(controlSchedules,...
        'RATE',10*meter^3/day,...
        'ORAT',10*meter^3/day,...
        'WRAT',10*meter^3/day,...
        'LRAT',10*meter^3/day,...
        'RESV',0,...
        'BHP',5*barsa));

    step = @(x0,u,varargin) mrstStep(x0,u,@mrstSimulationStep,stepSchedules(1),reservoirP,...
                                        'uScale',cellControlScales{1},...
                                        varargin{:});

    state = model.toStateVector( reservoirP.state);

    u  = schedules2CellControls( controlSchedules,'cellControlScales',cellControlScales);

    [ errorMax1 ] = testSimStepGradient(state,...
                                        u{1},...
                                        step,...
                                        'debug',true,...
                                        'pert',1e-5);

                                    
    ok = errorMax1 < 1e-4;

    mrstModule clear;
    warning('off','MATLAB:rmpath:DirNotFound');
    rmpath(genpath(remso_root));
    warning('on','MATLAB:rmpath:DirNotFound');
    
end

