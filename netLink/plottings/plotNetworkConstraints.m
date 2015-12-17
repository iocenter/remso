function [figN] =  plotNetworkConstraints(v, lbv, ubv, nScale, times, netCst, figN)
%%plotNetworkConstraints function that provide plots for network
%%constraints.

%% TODO: should receive size of each type of constraint, split the vector v, lbv, ubv, and call corresponding functions
    dp = zeros(netCst, numel(v));   
    ldp = zeros(netCst, numel(v));   
    udp  = zeros(netCst, numel(v));   
    
    for i=1:netCst
        figure(figN); 
        figN = figN+1;
        
        for j=1:numel(v)
            dp(i,j) = v{j}(end-netCst+i-1).*nScale(i);
            
            ldp(i,j) = lbv{j}(end-netCst+i-1).*nScale(i);            
            udp(i,j) = ubv{j}(end-netCst+i-1).*nScale(i);
        end
               
        
        plot(times.steps(2:end), dp(i,:), 'kx-', times.steps(2:end), ldp(i,:), 'rv-',times.steps(2:end), udp(i,:), 'r^-' );
        
        plot(times.steps(2:end), dp(i,:), 'kx-', times.steps(2:end), ldp(i,:), 'rv-',times.steps(2:end), udp(i,:), 'r^-' );

        title(strcat('Well: ', int2str(i)));       
        
        xlabel('time (days)');
        ylabel('frequency (Hertz)');
    end
end