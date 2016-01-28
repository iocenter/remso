function [figN] =  plotNetworkConstraints(v, lbv, ubv, nScale, times, netCst, figN, varargin)
%%plotNetworkConstraints function that provide plots for network
%%constraints.
    opt     = struct('freqCst', 0, 'pressureCst', 0, 'flowCst', 0, 'stepBreak', numel(v));
    opt     = merge_options(opt, varargin{:});


%% TODO: should receive size of each type of constraint, split the vector v, lbv, ubv, and call corresponding functions
%     dp = zeros(netCst, numel(v)); 
    
    dp = cell(numel(v), 1);
    ldp = cell(numel(v), 1);
    udp = cell(numel(v), 1);
    
    for i=1:numel(v)
%         if i<=opt.stepBreak
%             dp{i} = v{i}(end-opt.pressureCst:end-1).*nScale(end-opt.pressureCst+1:end);  %% TODO: scaling
%             ldp{i} = lbv{i}(end-opt.pressureCst:end-1).*nScale(end-opt.pressureCst+1:end);
%             udp{i} = ubv{i}(end-opt.pressureCst:end-1).*nScale(end-opt.pressureCst+1:end);
%         else
            dp{i} = v{i}(end-netCst:end-1).*nScale;
            ldp{i} = lbv{i}(end-netCst:end-1).*nScale;
            udp{i} = ubv{i}(end-netCst:end-1).*nScale;
%         end        
    end
    
%     for i=1:netCst
%         for j=1:numel(v{i})
%             dp(i,j) = v{j}(end-netCst+i-1).*nScale(i);
%             
%             ldp(i,j) = lbv{j}(end-netCst+i-1).*nScale(i);            
%             udp(i,j) = ubv{j}(end-netCst+i-1).*nScale(i);
%         end
%     end

    if opt.freqCst>0
        dpFreq  = zeros(opt.freqCst, numel(v));
        ldpFreq = zeros(opt.freqCst, numel(v));
        udpFreq = zeros(opt.freqCst, numel(v));       
%         for i=1:opt.stepBreak
%             dpFreq(:,i) = 0;
%             ldpFreq(:,i) = -inf;
%             udpFreq(:,i) = inf;
%         end
%         if opt.stepBreak<numel(v)            
            for j=1:numel(v)
                dpFreq(:,j) = dp{j}(1:opt.freqCst);
                ldpFreq(:,j) = ldp{j}(1:opt.freqCst);
                udpFreq(:,j) = udp{j}(1:opt.freqCst);
            end
%         end       
        
        figN = plotFrequencyConstraints(dpFreq, ldpFreq, udpFreq, times, opt.freqCst, figN);
    end

  
    if opt.pressureCst>0
        dpPressure  = zeros(opt.pressureCst, numel(v));
        ldpPressure = zeros(opt.pressureCst, numel(v));
        udpPressure = zeros(opt.pressureCst, numel(v));
        for i=1:numel(v)
            dpPressure(:,i) = dp{i}(end-opt.pressureCst+1:end);
            ldpPressure(:,i) = ldp{i}(end-opt.pressureCst+1:end);
            udpPressure(:,i) = udp{i}(end-opt.pressureCst+1:end);
        end
        
%         if opt.stepBreak<numel(v)
%             for j=opt.stepBreak+1:numel(v)
%                 dpPressure(:,j) = dp{j}(end-opt.freqCst+1:end);
%                 ldpPressure(:,j) = ldp{j}(end-opt.freqCst+1:end);
%                 udpPressure(:,j) = udp{j}(end-opt.freqCst+1:end);
%             end
%         end       
        
        figN = plotPressureConstraints(dpPressure, ldpPressure, udpPressure, times, opt.pressureCst, figN);
    end



    
    
    
    if opt.flowCst>0      
        figN = plotFlowConstraints(dp(opt.freqCst+opt.pressureCst+1:end, :), ldp(opt.freqCst+opt.pressureCst+1:end, :), udp(opt.freqCst+opt.pressureCst+1:end, :), times, opt.flowCst, figN);
    end
 
end