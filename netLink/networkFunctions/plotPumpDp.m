function [a ] = plotPumpDp(qf_start, qf_end, numFlows, freq_start, freq_end, numFreq, nStages, fref, qmin_60, qmax_60, varargin)
%plotPumpMap Plots the map of the ESPs (dp or dh versus flow rate)

    opt = struct('dp_map',true, 'dh_map',[]);
    opt = merge_options(opt, varargin{:});

    dh = zeros(numFlows,numFreq);    
    dp = zeros(numFlows, numFreq);
    flows = zeros(numFlows, numFreq);
    
    freq = zeros(numFreq, 1);
    interFlows = 1:numFlows;    
    interFreq  = 1:numFreq;
    for j=interFreq
        f = freq_start+ (j-1).*(freq_end-freq_start)./numFreq;
        freq(j)  =  f;   
        for i=interFlows
            q = qf_start+ (i-1).*(qf_end-qf_start)./numFlows;            

            flows(i,j) =  q./(meter^3/day);                    
            if opt.dp_map
                dp(i,j) = calcDp(q, f, fref, nStages);
            end
            if opt.dh_map
                dh(i,j)  = pumpDhExplicit(q)*((f./fref).^2)*nStages;
            end

        end    
    end

    a = figure(31);  
    
    colors = ['k-', 'r-', 'b-', 'm-', 'y-', 'kx-', 'rx-', 'bx-', 'mx-', 'yx-'];
    h = zeros(numFreq, 1);     
    for j=interFreq
        if opt.dp_map
            plot(flows, dp(:,j), colors(j));
        end
        
        if opt.dh_map
            plot(flows, dh(:,j), colors(j));
        end

%         %%TODO: display legend for each frequency    
%         legend(strcat('Freq: ', num2str(freq(j)), ' Hz,'));
        hold on;
    end
    
    title(strcat('Pump Map'));
    xlabel('flow (sm3/day)');    
    if opt.dp_map
        ylabel('Dp (bar)');    
    end
    if opt.dh_map
        ylabel('Dh (m)');      
    end
    
    qf_min = zeros(numFreq, 1);
    qf_max = zeros(numFreq, 1);        
    for j = interFreq(1:end-1)
        f1 = freq_start+ (j-1).*(freq_end-freq_start)./numFreq;
        f2 = freq_start+ (j).*(freq_end-freq_start)./numFreq;
        
        
        qf_min(j) = pump_rate(f1, qmin_60, 60);
        qf_min(j+1) = pump_rate(f2, qmin_60, 60);
        
        if opt.dh_map
            dp_min_1 = calcDh(qf_min(j),f1, fref, nStages);
            dp_min_2 = calcDh(qf_min(j+1),f2, fref, nStages);
        end
        
        if opt.dp_map        
            dp_min_1 = calcDp(qf_min(j),f1, fref, nStages);
            dp_min_2 = calcDp(qf_min(j+1),f2, fref, nStages);
        end
        
        line([qf_min(j)./(meter^3/day) qf_min(j+1)./(meter^3/day)], [dp_min_1,  dp_min_2]);
        
        qf_max(j) = pump_rate(f1, qmax_60, 60);
        qf_max(j+1) = pump_rate(f2, qmax_60, 60);
        
        if opt.dp_map
            dp_max_1 = calcDp(qf_max(j),f1, fref, nStages);
            dp_max_2 = calcDp(qf_max(j+1),f2, fref, nStages);
        end
        
        if opt.dh_map        
            dp_max_1 = calcDh(qf_max(j),f1, fref, nStages);
            dp_max_2 = calcDh(qf_max(j+1),f2, fref, nStages);
        end
        
        line([qf_max(j)./(meter^3/day) qf_max(j+1)./(meter^3/day)], [dp_max_1,  dp_max_2]);
       
    end
end


function dh = calcDh(q,f,fref,nStages)
    dh  = pumpDhExplicit(q)*((f./fref).^2)*nStages;
end

function dp = calcDp(q, f, fref, nStages)        
    dh  = pumpDhExplicit(q)*((f./fref).^2)*nStages;
    mixDens = 900;
    dp  = (dh*norm(gravity)*mixDens)./barsa;

end
