function [ ] = plotPumpMap(netSol, wellSol, forwardStates, schedule, p, nScale, numStages, pScale, wellId)
%plotPumpMap Plots the map of the ESPs (dh versus flow rate)
   
    qf_start = 5*(meter^3/day);
    qf_end   = 500.*(meter^3/day);
    numIntervals = 100;
  
    dh = zeros(numIntervals,1);
    freq = zeros(numIntervals,1);
    flows = zeros(numIntervals,1);
    interval = 1:numIntervals;    
    for i=interval
        q = qf_start+ (i-1).*(qf_end-qf_start)./numIntervals;
        
        flows(i) =  q;
        
        [dh(i) freq(i)]  = pumpDh(netSol, wellSol, forwardStates, p, numStages, pScale, wellId, q);        
        
    end    
%      plot(times.steps(2:end), dh(i,:), 'kx-', times.steps(2:end), ldp(i,:), 'rv-',times.steps(2:end), udp(i,:), 'r^-' );
%      plot(times.steps(2:end), dh(i,:), 'kx-', times.steps(2:end), ldp(i,:), 'rv-',times.steps(2:end), udp(i,:), 'r^-' );

    figure(31); 
      title(strcat('Pump Map: Well', int2str(wellId)));
     plot(flows./(meter^3/day), dh, 'kx-');
     xlabel('flow (sm3/day)');
     ylabel('Dh (m)');
     
     figure(32); 
     
     plot(flows./(meter^3/day), freq,'kx-');
        
     title(strcat('Pump Frequency: Well', int2str(wellId)));
     
     xlabel('flow (sm3/day)');
     ylabel('Frequency (Hz)');
    
     
%      qf_min = pump_rate(p, qf.*0, 60);
%      qf_max   = pump_rate(p, qf.*10, 60);

end


function [dhf freq] = pumpDh(netSol, wellSol, forwardStates,   p, numStages, pScale, wellId, q)
    vw = getVertex(netSol, wellId);        
    ew = getEdge(netSol, vertcat(vw.Eout));

    wcut = vw.qwV./(vw.qoV+vw.qwV);
        
    vw.qoV = -q.*(1-wcut);
    vw.qwV = -q.*(wcut);
        
    netSol = updateVertex(netSol, vw);
    netSol = runNetwork(netSol, wellSol, forwardStates{numel(forwardStates)}, p, pScale, 'sensitivityAnalysis', true);   % running the network
    
   
    dpf = getChokesDp(netSol); % dp in the pumps    
    dpPump = dpf(1);    
   
    inletStr = vertcat(ew.stream);
    mixtureDen = (vertcat(inletStr.oil_dens) + vertcat(inletStr.water_dens))./2;  % density of the mixture
        
    dhf= pump_dh(dpPump, mixtureDen); % dh in the pumps
    
    freq = pump_eq_system_explicit(q, dhf, 60, numStages);  % solves a system of equations to obtain frequency, flow and dh at 60Hz
    
    if isa(freq, 'ADI')
        freq.val = real(freq.val);
        freq.jac = cellfun(@(w) real(w) ,freq.jac, 'UniformOutput', false);
    else
        freq = real(freq);
    end
    
    freq = freq(1);
    
end



