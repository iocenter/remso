function [ ] = plotPumpDh(qf_start, qf_end, numFlows, freq_start, freq_end, numFreq, nStages, fref)
%plotPumpMap Plots the map of the ESPs (dh versus flow rate)
    dh = zeros(numFlows,numFreq);    
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

        end    
    end
%      plot(times.steps(2:end), dh(i,:), 'kx-', times.steps(2:end), ldp(i,:), 'rv-',times.steps(2:end), udp(i,:), 'r^-' );
%      plot(times.steps(2:end), dh(i,:), 'kx-', times.steps(2:end), ldp(i,:), 'rv-',times.steps(2:end), udp(i,:), 'r^-' );

    figure(31);
    title(strcat('Pump Map'));
    xlabel('flow (sm3/day)');
    ylabel('Dh (m)');
    colors = ['k-', 'r-', 'b-', 'm-', 'y-', 'kx-', 'rx-', 'bx-', 'mx-', 'yx-'];
    h = zeros(numFreq, 1); 
    
    
    for j=interFreq
        plot(flows./(meter^3/day), dh(:,j), colors(j));           
        
        legend(strcat('Freq: ', num2str(freq(j)), ' Hz,'));
        hold on;
    end
    
    title(strcat('Pump Map'));
    xlabel('flow (sm3/day)');
    ylabel('Dh (m)');      
     
    
    qf_min = zeros(numFreq, 1);
    qf_max = zeros(numFreq, 1);
    eps = 10;
    
    for j = interFreq
       qf_min(j) = pump_rate(freq(j), 636, 60);
       qf_max(j) = pump_rate(freq(j), 1542, 60);       
       
       if j<numFreq
            line([qf_min(j) qf_min(j)],[dh(1,j) dh(1,j)+eps]);           
            line([qf_max(j) qf_max(j)],[dh(1,j) dh(1,j)+eps]);
       end
    end
end