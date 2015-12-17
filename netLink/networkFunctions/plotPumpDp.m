function [ ] = plotPumpDp(qf_start, qf_end, numFlows, freq_start, freq_end, numFreq, nStages, fref)
%plotPumpMap Plots the map of the ESPs (dp versus flow rate)
    dh = zeros(numFlows,numFreq);    
    dp = zeros(numFlows, numFreq);
    flows = zeros(numFlows, numFreq);
    
    freq = zeros(numFreq, 1);
    interFlows = 1:numFlows;    
    interFreq  = 1:numFreq;
    for j=interFreq
        f = freq_start+ (j).*(freq_end-freq_start)./numFreq;
        freq(j)  =  f;   
        for i=interFlows
            q = qf_start+ (i).*(qf_end-qf_start)./numFlows;            

            flows(i,j) =  q;
                     

            dh(i,j)  = pumpDhExplicit(q)*((f./fref).^2)*nStages;
            mixDens = 900;
            dp(i,j)  = (dh(i,j)*norm(gravity)*mixDens)./barsa;

        end    
    end

    figure(31);  
    
    colors = ['k-', 'r-', 'b-', 'm-', 'y-', 'kx-', 'rx-', 'bx-', 'mx-', 'yx-'];
    h = zeros(numFreq, 1);     
    for j=interFreq
        plot(flows./(meter^3/day), dp(:,j), colors(j));           

%         %%TODO: display legend for each frequency    
%         legend(strcat('Freq: ', num2str(freq(j)), ' Hz,'));
        hold on;
    end
    
    title(strcat('Pump Map'));
    xlabel('flow (sm3/day)');
    ylabel('Dp (m)');      
     
    
    qf_min = zeros(numFreq, 1);
    qf_max = zeros(numFreq, 1);
    eps = 1;
    
    for j = interFreq
%        qf_min(j) = pump_rate(freq(j), 636, 60);
%        qf_max(j) = pump_rate(freq(j), 1542, 60);       
        
       qf_min(j) = pump_rate(freq(j), 5, 60);
       qf_max(j) = pump_rate(freq(j), 500, 60);       
       
       if j<numFreq
            line([qf_min(j) qf_min(j)],[dp(1,j) dp(1,j)+eps]);           
            line([qf_max(j) qf_max(j)],[dp(1,j) dp(1,j)+eps]);
       end
    end
end