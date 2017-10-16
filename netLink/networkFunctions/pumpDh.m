function [ dh ] = pumpDh(q)
%calculates a single stage dh of the pump for the flow rate q at 60 Hz 

    a0 = 19.37;       % m 
    a1 = -1.23e-02;     % in m/sm3/d
    a2 = 2.24e-05;      % in m/(sm3/d)^2      
    a3 = -1.86e-08;     % in m/(sm3/d)^3
    a4 = 4.13e-12;      % in m/(sm3/d)^4  
    
    dh = a4.*((q.^4)./(meter^3/day)) +  a3.*((q.^3)./(meter^3/day)) +  a2.*((q.^2)./(meter^3/day)) +  a1.*(q./(meter^3/day)) + a0;
    
end

