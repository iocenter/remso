function [figN] =  plotPressureConstraints(dp, ldp, udp, times, netCst, figN)
%%plotFrequencyConstraints plots constraints on the pressure difference in the pumps
%%(ESPs) of the network    
    for i=1:netCst
        figure(figN); 
        figN = figN+1;
        
        plot(times.steps(2:end), dp(i,:)./barsa, 'kx-', times.steps(2:end), ldp(i,:)./barsa, 'rv-',times.steps(2:end), udp(i,:)./barsa, 'r^-' );
        
        plot(times.steps(2:end), dp(i,:)./barsa, 'kx-', times.steps(2:end), ldp(i,:)./barsa, 'rv-',times.steps(2:end), udp(i,:)./barsa, 'r^-' );

        title(strcat('Well: ', int2str(i)));       
        
        xlabel('time (days)');
        ylabel('dp (bar)');
    end
end