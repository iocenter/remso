function [figN] =  plotNetworkConstraints(v, lbv, ubv, nScale, times, netCst, figN, varargin)
%%plotNetworkConstraints function that provide plots for network
%%constraints.
    opt     = struct('freqCst', 0, ...
                    'pressureCst', 0, ...
                    'flowCst', 0, ...
                    'stepBreak', numel(v), ...
                    'nW', 7, ... 
                    'vScale', [], ...
                    'extremePoints', [], ...
                    'qlMin', [], ...
                    'qlMax', [], ...
                    'nStages', [], ...
                    'freqMin', [], ...
                    'freqMax', [], ...
                    'baseFreq', []);
    opt     = merge_options(opt, varargin{:});

    n = cellfun(@(vk) vk(end-netCst:end-1), v, 'UniformOutput', false); %%     
    lbn  = cellfun(@(vk) vk(end-netCst:end-1), lbv, 'UniformOutput', false); %%     
    ubn  = cellfun(@(vk) vk(end-netCst:end-1), ubv, 'UniformOutput', false); %%     
    
    qf = cell(numel(v), 1);
    wcut = cell(numel(v),1);
    dp = cell(numel(v), 1);
    ldp = cell(numel(v), 1);
    udp = cell(numel(v), 1);
    
    for i=1:numel(v)
        qf{i} = v{i}(3:opt.nW).*opt.vScale(3:opt.nW) + v{i}(opt.nW+3:2*opt.nW).*opt.vScale(opt.nW+3:2*opt.nW);  %% liquid flow rate
        wcut{i} = (v{i}(3:opt.nW).*opt.vScale(3:opt.nW))./qf{i};
        
        dp{i}  = n{i}.*nScale;
        ldp{i} = lbn{i}.*nScale;
        udp{i} = ubn{i}.*nScale;

    end    
    if ~isempty(opt.extremePoints)
        qfFlows = zeros(opt.pressureCst, numel(v));
        dpPressure  = zeros(opt.pressureCst, numel(v));
        ldpPressure = zeros(opt.pressureCst, numel(v));
        udpPressure = zeros(opt.pressureCst, numel(v));
        for j=1:numel(v)
            qfFlows(:,j) = -qf{j};
            dpPressure(:,j) = dp{j}(end-opt.pressureCst+1:end);
            ldpPressure(:,j) = ldp{j}(end-opt.pressureCst+1:end);
            udpPressure(:,j) = udp{j}(end-opt.pressureCst+1:end);
        end
        figN = plotLinearPumpConstraints(qfFlows, dpPressure, times, opt.pressureCst, figN, 'extremePoints', opt.extremePoints, 'qlMin', opt.qlMin, 'qlMax', opt.qlMax, 'nStages', opt.nStages);        
    else
        if opt.flowCst>0
            dpFlows = zeros(opt.flowCst/2, numel(v));
            wcutFlows = zeros(opt.flowCst/2, numel(v));
            ldpFlows = zeros(opt.flowCst/2, numel(v));
            udpFlows = zeros(opt.flowCst/2, numel(v));
            
            
            for i=1:numel(v)
                wcutFlows(:,i) = wcut{i};
                dpFlows(:,i) = -qf{i};
                minFlow = dp{i}(1:opt.flowCst/2);  %% related to the constraints called minFlow
                maxFlow = dp{i}(opt.flowCst/2+1:opt.flowCst);   %% related to the constraints called maxFlow
                
                ldpFlows(:,i) = -minFlow - qf{i};   %% this calculation gives the min flow rate at the desired frequency
                udpFlows(:,i) = maxFlow - qf{i};    %% this calculation gives the max flow rate at the desired frequency
            end
        end        
        
        if opt.freqCst>0
            dpFreq  = zeros(opt.freqCst, numel(v));
            ldpFreq = zeros(opt.freqCst, numel(v));
            udpFreq = zeros(opt.freqCst, numel(v));
            
            for j=1:numel(v)
                dpFreq(:,j) = dp{j}(opt.flowCst+1:opt.flowCst+opt.freqCst);
                ldpFreq(:,j) = ldp{j}(opt.flowCst+1:opt.flowCst+opt.freqCst);
                udpFreq(:,j) = udp{j}(opt.flowCst+1:opt.flowCst+opt.freqCst);
            end
            
        end         
        figN = plotFrequencyConstraints(dpFreq, ldpFreq, udpFreq, times, opt.freqCst, figN);
        
        if opt.pressureCst>0
            dpPressure  = zeros(opt.pressureCst, numel(v));
            ldpPressure = zeros(opt.pressureCst, numel(v));
            udpPressure = zeros(opt.pressureCst, numel(v));
            for i=1:numel(v)
                dpPressure(:,i) = dp{i}(end-opt.pressureCst+1:end);
                ldpPressure(:,i) = ldp{i}(end-opt.pressureCst+1:end);
                udpPressure(:,i) = udp{i}(end-opt.pressureCst+1:end);
            end           
        end   

        % bounds for pump frequencies
        qminFmin = pump_rate(opt.freqMin, opt.qlMin, opt.baseFreq);
        qminFmax = pump_rate(opt.freqMax, opt.qlMin, opt.baseFreq);

        qmaxFmin = pump_rate(opt.freqMin, opt.qlMax, opt.baseFreq);
        qmaxFmax = pump_rate(opt.freqMax, opt.qlMax, opt.baseFreq);

        qf_start = [qminFmin, qminFmax];
        qf_end = [qmaxFmin, qmaxFmax];        
        for i=1:numel(opt.nStages)
            figure(figN);
            figN = figN + 1;
            
            plotNonlinearPumpMap(dpFlows(i,:), dpPressure(i,:), times, qf_start(i,:), qf_end(i,:), 10, opt.freqMin(i,:), opt.freqMax(i,:), 2, opt.nStages(i,:), opt.baseFreq(i,:), opt.qlMin(i,:), opt.qlMax(i,:));   
            
            title(strcat('Pump Map: p ', int2str(i)));
            xlabel('Flow rate (sm3/day)');
            ylabel('Pressure Difference (bar)');       
        end        
    end
end