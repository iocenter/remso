function [figN] =  plotFrequencyConstraints(freq, lfreq, ufreq, times, netCst, figN)
%%plotFrequencyConstraints plots constraints on the frequency of pumps
%%(ESPs) of the network    
    for i=1:netCst
        figure(figN); 
        figN = figN+1;
        
        plot(times.steps(2:end), freq(i,:), 'kx-', times.steps(2:end), lfreq(i,:), 'rv-',times.steps(2:end), ufreq(i,:), 'r^-' );
        
        plot(times.steps(2:end), freq(i,:), 'kx-', times.steps(2:end), lfreq(i,:), 'rv-',times.steps(2:end), ufreq(i,:), 'r^-' );

        title(strcat('Well: ', int2str(i)));       
        
        xlabel('time (days)');
        ylabel('frequency (Hertz)');
    end
end