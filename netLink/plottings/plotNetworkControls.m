function [figN] =  plotNetworkControls(u, lbu, ubu, uScale, times, numControls, figN)
   
    pSep = zeros(numControls, numel(u));   
    ldp = zeros(numControls, numel(u));   
    udp  = zeros(numControls, numel(u));   
    
    for i=1:numControls
        figure(figN); 
        figN = figN+1;
        
        for j=1:numel(u)
            
            pSep(i,j) = u{j}(end-numControls+i).*uScale{j}(end-numControls+i);   % frequency in Hz
            
            ldp(i,j) = lbu{j}(end-numControls+i).*uScale{j}(end-numControls+i);            
            udp(i,j) = ubu{j}(end-numControls+i).*uScale{j}(end-numControls+i);
        end
               
        
        plot(times.controls(2:end), pSep(i,:), 'kx-', times.controls(2:end), ldp(i,:), 'rv-',times.controls(2:end), udp(i,:), 'r^-' );

        title(strcat('Frequency (Hz)'));       
        
        xlabel('time (day)');
        ylabel('freq (Hz)');
    end
end

