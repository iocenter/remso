clear; clc;

here = fileparts(mfilename('fullpath'));
if isempty(here)
here = pwd();
end

addpath(genpath(fullfile(here,filesep,'..',filesep,'..',filesep,'mrstDerivated')));
addpath(genpath(fullfile(here,filesep,'..',filesep,'..',filesep,'mrstLink')));
addpath(genpath(fullfile(here,filesep,'..',filesep,'..',filesep,'mrstLink',filesep,'wrappers',filesep,'procedural')));
addpath(genpath(fullfile(here,filesep,'..',filesep,'..',filesep,'netLink')));
addpath(genpath(fullfile(here,filesep,'..',filesep,'..',filesep,'netLink',filesep,'plottings')));
addpath(genpath(fullfile(here,filesep,'..',filesep,'..',filesep,'optimization',filesep,'multipleShooting')));
addpath(genpath(fullfile(here,filesep,'..',filesep,'..',filesep,'optimization',filesep,'plotUtils')));
addpath(genpath(fullfile(here,filesep,'..',filesep,'..',filesep,'optimization',filesep,'remso')));
addpath(genpath(fullfile(here,filesep,'..',filesep,'..',filesep,'optimization',filesep,'remsoSequential')));
addpath(genpath(fullfile(here,filesep,'..',filesep,'..',filesep,'optimization',filesep,'remsoCrossSequential')));
addpath(genpath(fullfile(here,filesep,'..',filesep,'..',filesep,'optimization',filesep,'singleShooting')));
addpath(genpath(fullfile(here,filesep,'..',filesep,'..',filesep,'optimization',filesep,'utils')));
addpath(genpath(fullfile(here,filesep,'reservoirData')));

% qf_start = 636*(meter^3/day);
% qf_end   = 1542.*(meter^3/day);

qf_start = 15*(meter^3/day);
qf_end   = 900.*(meter^3/day);
numFlows = 50;  
    
freq_start = 20;
freq_end  = 90;
numFreq  = 10;

numStages = 80;
fref = 60;

qmin_60 = 30*(meter^3/day);
qmax_60 = 600*(meter^3/day);

% plotPumpDh(qf_start, qf_end, numFlows, freq_start, freq_end, numFreq, numStages, fref);

a = plotPumpDp(qf_start, qf_end, numFlows, freq_start, freq_end, numFreq, numStages, fref, qmin_60, qmax_60, 'dp_map', true, 'dh_map', false);

% load optFiles/greedyStrategy.mat;
% flowScale = 10*meter^3/day;
% pressureScale = 5*barsa;
% for i=1:numel(v)
%     ql = flowScale*(abs((v{i}(7) + v{i}(7+7))))/(meter^3/day);
%     
%     
% %     dp = pressureScale*(v{i}(end-5))./barsa;
%     dp = abs(min(0,pressureScale*(v{i}(end-5))./barsa));    
%     plot(ql, dp, '*r');
%     
%     saveas(a, strcat('p5/',num2str(i),'.png'));
% end


