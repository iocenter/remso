function [figN] =  plotNetworkConstraints(v, lbv, ubv, nScale, times, netCst, figN, varargin)
%%plotNetworkConstraints function that provide plots for network
%%constraints.
    opt     = struct('freqCst', 0, 'pressureCst', 0, 'flowCst', 0, 'stepBreak', numel(v));
    opt     = merge_options(opt, varargin{:});
    
    dp = cell(numel(v), 1);
    ldp = cell(numel(v), 1);
    udp = cell(numel(v), 1);
    
    for i=1:numel(v)
        if i<=opt.stepBreak
            dp{i} = v{i}(end-netCst:end-1).*nScale;
            ldp{i} = lbv{i}(end-netCst:end-1).*nScale;
            udp{i} = ubv{i}(end-netCst:end-1).*nScale;
        else
            dp{i} = v{i}(end-netCst:end-1).*nScale;
            ldp{i} = lbv{i}(end-netCst:end-1).*nScale;
            udp{i} = ubv{i}(end-netCst:end-1).*nScale;
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Plot flow constraints %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if opt.flowCst>0
        dpFlow  = zeros(opt.flowCst, numel(v));
        ldpFlow = zeros(opt.flowCst, numel(v));
        udpFlow = zeros(opt.flowCst, numel(v));
        
        for i=1:numel(v)
            dpFlow(:,i) = dp{i}(1:opt.flowCst);
            ldpFlow(:,i) = ldp{i}(1:opt.flowCst);
            udpFlow(:,i) = udp{i}(1:opt.flowCst);
        end
        
        figN = plotFlowConstraints(dpFlow, ldpFlow, udpFlow, times, opt.flowCst, figN);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Plot Frequency constraints %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if opt.freqCst>0
        dpFreq  = zeros(opt.freqCst, numel(v));
        ldpFreq = zeros(opt.freqCst, numel(v));
        udpFreq = zeros(opt.freqCst, numel(v));
        %         for i=1:opt.stepBreak
        %             dpFreq(:,i) = dp{i}(1:opt.freqCst);
        %             ldpFreq(:,i) = lbv{i}(1:opt.freqCst);
        %             udpFreq(:,i) = ubv{i}(1:opt.freqCst);
        %         end
        %         if opt.stepBreak<numel(v)
        for j=1:numel(v)
            dpFreq(:,j) = dp{j}(opt.flowCst+1:opt.flowCst+opt.freqCst);
            ldpFreq(:,j) = ldp{j}(opt.flowCst+1:opt.flowCst+opt.freqCst);
            udpFreq(:,j) = udp{j}(opt.flowCst+1:opt.flowCst+opt.freqCst);
        end        
        figN = plotFrequencyConstraints(dpFreq, ldpFreq, udpFreq, times, opt.freqCst, figN);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Plot pressure constraints %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if opt.pressureCst>0
        dpPressure  = zeros(opt.pressureCst, numel(v));
        ldpPressure = zeros(opt.pressureCst, numel(v));
        udpPressure = zeros(opt.pressureCst, numel(v));
        for i=1:numel(v)
            dpPressure(:,i) = dp{i}(end-opt.pressureCst+1:end);
            ldpPressure(:,i) = ldp{i}(end-opt.pressureCst+1:end);
            udpPressure(:,i) = udp{i}(end-opt.pressureCst+1:end);
        end
        
        figN = plotPressureConstraints(dpPressure, ldpPressure, udpPressure, times, opt.pressureCst, figN);
    end
    

   
end