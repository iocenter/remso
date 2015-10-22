function [figN] =  plotNetworkConstraints(v, lbv, ubv, nScale, times, netCst, figN)
   
    dp = zeros(netCst, numel(v));   
    ldp = zeros(netCst, numel(v));   
    udp  = zeros(netCst, numel(v));   
    
    for i=1:netCst
        figure(figN); 
        figN = figN+1;
        
        for j=1:numel(v)
            dp(i,j) = v{j}(end-netCst+i-1).*nScale(i)/barsa; % dp in equipment (excludes last element since it is the NPV)
            
            ldp(i,j) = lbv{j}(end-netCst+i-1).*nScale(i)/barsa;            
            udp(i,j) = ubv{j}(end-netCst+i-1).*nScale(i)/barsa;
        end
               
        
        plot(times.steps(2:end), dp(i,:), 'kx-', times.steps(2:end), ldp(i,:), 'rv-',times.steps(2:end), udp(i,:), 'r^-' );
        
%         plot(times.steps(2:end), ldp(i,:), 'rv-');
        
        plot(times.steps(2:end), dp(i,:), 'kx-', times.steps(2:end), ldp(i,:), 'rv-',times.steps(2:end), udp(i,:), 'r^-' );
        
%         plot(times.steps(2:end), ldp(i,:), 'rv-');
        
%         plot(times.steps(2:end), udp(i,:), 'r^-');

        title(strcat('Well: ', int2str(i)));       
        
        xlabel('time (sec)');
        ylabel('dp in equipment (bar)');
    end
end

