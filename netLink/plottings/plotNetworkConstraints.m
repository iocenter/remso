function [figN] =  plotNetworkConstraints(v, lbv, ubv, nScale, times, netCst, figN, varargin)
%%plotNetworkConstraints function that provide plots for network
%%constraints.
    opt     = struct('freqCst', 0, 'pressureCst', 0, 'flowCst', 0, 'stepBreak', numel(v), 'nW', 7, 'vScale', []);
    opt     = merge_options(opt, varargin{:});


%% TODO: should receive size of each type of constraint, split the vector v, lbv, ubv, and call corresponding functions
%     dp = zeros(netCst, numel(v)); 

    n = cellfun(@(vk) vk(end-netCst:end-1), v, 'UniformOutput', false); %%     
    lbn  = cellfun(@(vk) vk(end-netCst:end-1), lbv, 'UniformOutput', false); %%     
    ubn  = cellfun(@(vk) vk(end-netCst:end-1), ubv, 'UniformOutput', false); %%     
    
    qf = cell(numel(v), 1);
    dp = cell(numel(v), 1);
    ldp = cell(numel(v), 1);
    udp = cell(numel(v), 1);
    
    for i=1:numel(v)
%         if i<=opt.stepBreak
%             dp{i} = v{i}(end-opt.pressureCst:end-1).*nScale(end-opt.pressureCst+1:end);  %% TODO: scaling
%             ldp{i} = lbv{i}(end-opt.pressureCst:end-1).*nScale(end-opt.pressureCst+1:end);
%             udp{i} = ubv{i}(end-opt.pressureCst:end-1).*nScale(end-opt.pressureCst+1:end);
%         else
            qf{i} = v{i}(3:opt.nW).*opt.vScale(3:opt.nW) + v{i}(opt.nW+3:2*opt.nW).*opt.vScale(opt.nW+3:2*opt.nW);  %% liquid flow rate             
            
            dp{i}  = n{i}.*nScale; 
            ldp{i} = lbn{i}.*nScale;
            udp{i} = ubn{i}.*nScale;
%         end        
    end
    
    %% TODO: plot a single graphic for flow constraints
    if opt.flowCst>0   
        dpFlows = zeros(opt.flowCst/2, numel(v));
        ldpFlows = zeros(opt.flowCst/2, numel(v));
        udpFlows = zeros(opt.flowCst/2, numel(v));
        
        
        for i=1:numel(v)
            dpFlows(:,i) = -qf{i};
            minFlow = dp{i}(1:opt.flowCst/2);  %% related to the constraints called minFlow                        
            maxFlow = dp{i}(opt.flowCst/2+1:opt.flowCst);   %% related to the constraints called maxFlow 
            
            ldpFlows(:,i) = -minFlow - qf{i};   %% this calculation gives the min flow rate at the desired frequency            
            udpFlows(:,i) = maxFlow - qf{i};    %% this calculation gives the max flow rate at the desired frequency  
        end
        
        figN = plotFlowConstraints(dpFlows, ldpFlows, udpFlows, times, opt.flowCst./2, figN);
    end

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
                dpFreq(:,j) = dp{j}(opt.flowCst+1:opt.flowCst+opt.freqCst);
                ldpFreq(:,j) = ldp{j}(opt.flowCst+1:opt.flowCst+opt.freqCst);
                udpFreq(:,j) = udp{j}(opt.flowCst+1:opt.flowCst+opt.freqCst);
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
    
    
  
 
end