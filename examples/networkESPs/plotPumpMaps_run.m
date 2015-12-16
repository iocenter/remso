% qf_start = 636*(meter^3/day);
% qf_end   = 1542.*(meter^3/day);

qf_start = 300*(meter^3/day);
qf_end   = 1842.*(meter^3/day);
numFlows = 50;  

freq_start = 30;
freq_end   = 60;
numFreq  = 10;

numStages = 100;
fref = 60;

plotPumpDh(qf_start, qf_end, numFlows, freq_start, freq_end, numFreq, numStages, fref);
    
