mrstModule clear
mrstModule add deckformat ad-core ad-blackoil ad-props co2lab optimization

here = fileparts(mfilename('fullpath'));
if isempty(here)
    here = pwd();
end
remso_root = fullfile(here,filesep,'..',filesep,'..',filesep);

addpath(genpath(fullfile(remso_root,'mrstDerivated')));
addpath(genpath(fullfile(remso_root,'mrstLink',filesep,'utils')));
addpath(genpath(fullfile(remso_root,'mrstLink',filesep,'wrappers',filesep,'procedural')));
addpath(genpath(fullfile(remso_root,'mrstLink',filesep,'wrappers',filesep,'OOP')));
addpath(genpath(fullfile(remso_root,'optimization',filesep,'testFunctions')));

addpath(genpath(fullfile(here,filesep,'reservoirData')));

[reservoirP] = initReservoir( 'reallySimpleRes.data');

model = reservoirP.model;
state = reservoirP.state;
schedule = reservoirP.schedule;

% make sure the tolerances are tight!
model.toleranceMB = min(1e-12,model.toleranceMB);
model.toleranceCNV = min(1e-7,model.toleranceCNV);
model.toleranceWellBHP = min(barsa/1e-5,model.toleranceWellBHP);
model.toleranceWellRate = min(1/day/1e-5,model.toleranceWellRate);    

NPVScale = 1e3;
target = @(wellSols,states,varargin) NPVOW(model.G,...
                                           wellSols,...
                                           states,...
                                           'OilPrice',             1.0/NPVScale, ...
                                           'WaterProductionCost',  0.1/NPVScale, ...
                                           'WaterInjectionCost',   0.1/NPVScale,...
                                           varargin{:});

ok = testTargetGradient(state,model,schedule,target);


%warning('off','MATLAB:rmpath:DirNotFound');
%mrstModule clear
%rmpath(genpath(fullfile(remso_root,'mrstDerivated')));
%rmpath(genpath(fullfile(remso_root,'mrstLink',filesep,'utils')));
%rmpath(genpath(fullfile(remso_root,'mrstLink',filesep,'wrappers',filesep,'procedural')));
%rmpath(genpath(fullfile(remso_root,'mrstLink',filesep,'wrappers',filesep,'OOP')));
%rmpath(genpath(fullfile(remso_root,'optimization',filesep,'testFunctions')));
%rmpath(genpath(fullfile(here,filesep,'reservoirData')));
%warning('on','MATLAB:rmpath:DirNotFound');
