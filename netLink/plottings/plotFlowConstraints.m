function [figN] =  plotFlowConstraints(qf, lqf, uqf, times, netCst, figN)
%%plotFrequencyConstraints plots constraints on the flows in the pumps
%%(ESPs) of the network    
    for i=1:netCst
        figure(figN); 
        figN = figN+1;
        
        plot(times.steps(2:end), qf(i,:)./(meter^3/day), 'kx-', times.steps(2:end), lqf(i,:)./(meter^3/day), 'rv-',times.steps(2:end), uqf(i,:)./(meter^3/day), 'r^-' );
        
        plot(times.steps(2:end), qf(i,:)./(meter^3/day), 'kx-', times.steps(2:end), lqf(i,:)./(meter^3/day), 'rv-',times.steps(2:end), uqf(i,:)./(meter^3/day), 'r^-' );

        title(strcat('Well: ', int2str(i)));       
        
        xlabel('time (days)');
        ylabel('flow rate (sm3/d)');
    end
end