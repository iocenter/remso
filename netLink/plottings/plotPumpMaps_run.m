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

qf_start = 2.5*(meter^3/day);
qf_end   = 350.*(meter^3/day);  
numFlows = 50;  
    
freq_start = 30;
freq_end  = 90;
numFreq  = 10;

numStages = 100;
fref = 60;

qmin_60 = 5*(meter^3/day);
qmax_60 = 180*(meter^3/day);

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

load itVars.mat
numInjectors = 2;
numProducers = 5;
prodIndex = 4;

flowScale = 10*meter^3/day;
pressureScale = 5*barsa;
 for i=1:numel(v)     
     ql = flowScale*(abs((v{i}(numInjectors+prodIndex) + v{i}((numInjectors+numProducers)+numInjectors+prodIndex))))/(meter^3/day);  
     dp = abs(min(0,pressureScale*(v{i}(end-numProducers))./barsa));    
%      if ql > 100
%          i  
%      end
     if i<16
        plot(ql, dp, '*r');    
     end
     
 end




