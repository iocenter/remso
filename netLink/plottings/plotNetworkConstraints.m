function [figN] =  plotNetworkConstraints(v, lbv, ubv, nScale, times, netCst, figN, varargin)
%%plotNetworkConstraints function that provide plots for network
%%constraints.
    opt     = struct('freqCst', 0, 'pressureCst', 0, 'flowCst', 0);
    opt     = merge_options(opt, varargin{:});


%% TODO: should receive size of each type of constraint, split the vector v, lbv, ubv, and call corresponding functions
    dp = zeros(netCst, numel(v));   
    ldp = zeros(netCst, numel(v));   
    udp  = zeros(netCst, numel(v));   
    
    
    for i=1:netCst
        for j=1:numel(v)
            dp(i,j) = v{j}(end-netCst+i-1).*nScale(i);
            
            ldp(i,j) = lbv{j}(end-netCst+i-1).*nScale(i);            
            udp(i,j) = ubv{j}(end-netCst+i-1).*nScale(i);
        end
    end
            
    if opt.freqCst>0
        figN = plotFrequencyConstraints(dp(1:opt.freqCst, :), ldp(1:opt.freqCst, :), udp(1:opt.freqCst, :), times, opt.freqCst, figN);
    end
    
    if opt.pressureCst>0
        figN = plotPressureConstraints(dp(opt.freqCst+1:opt.freqCst+opt.pressureCst, :), ldp(opt.freqCst+1:opt.freqCst+opt.pressureCst, :), udp(opt.freqCst+1:opt.freqCst+opt.pressureCst, :), times, opt.pressureCst, figN);
    end
    
    if opt.flowCst>0      
        figN = plotFlowConstraints(dp(opt.freqCst+opt.pressureCst+1:end, :), ldp(opt.freqCst+opt.pressureCst+1:end, :), udp(opt.freqCst+opt.pressureCst+1:end, :), times, opt.flowCst, figN);
    end
 
end