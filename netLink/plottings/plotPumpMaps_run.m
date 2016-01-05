clear; clc;

addpath(genpath('../../mrstDerivated'));
addpath(genpath('../../mrstLink'));
addpath(genpath('../../mrstLink/wrappers/procedural'));
addpath(genpath('../../netLink'));
addpath(genpath('../../netLink/plottings'));
addpath(genpath('../../optimization/multipleShooting'));
addpath(genpath('../../optimization/plotUtils'));
addpath(genpath('../../optimization/remso'));
addpath(genpath('../../optimization/remsoSequential'));
addpath(genpath('../../optimization/remsoCrossSequential'));
addpath(genpath('../../optimization/singleShooting'));
addpath(genpath('../../optimization/utils'));
addpath(genpath('reservoirData'));

% qf_start = 636*(meter^3/day);
% qf_end   = 1542.*(meter^3/day);

qf_start = 1*(meter^3/day);
qf_end   = 900.*(meter^3/day);
numFlows = 50;  

freq_start = 30;
freq_end  = 60;
numFreq  = 20;

numStages = 50;
fref = 60;

qmin_60 = 50*(meter^3/day);
qmax_60 = 600*(meter^3/day);

% plotPumpDh(qf_start, qf_end, numFlows, freq_start, freq_end, numFreq, numStages, fref);

plotPumpDp(qf_start, qf_end, numFlows, freq_start, freq_end, numFreq, numStages, fref, qmin_60, qmax_60);
    
