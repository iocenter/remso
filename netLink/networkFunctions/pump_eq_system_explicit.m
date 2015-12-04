function [freq ] = pump_eq_system_explicit(qf, dhf, fref, nStages)
% Calculates the head difference for the pump at a given pressure loss.
% It can be used to compute upper and lower bounds for the pressure loss
% in the equipment. The coefficients of the curve were found in the paper
% 'Exploring the potential of model-based optimization in oil production
% gathering networks with esp-produced, high water cut wells'.
% A system of equations is solved to obtain the frequency of the pump for a 
% given pressure drop and flow rate flowing through it.
    a0 = 19.37;       % m 
    a1 = -1.23e-02;     % in m/sm3/d
    a2 = 2.24e-05;      % in m/(sm3/d)^2      
    a3 = -1.86e-08;     % in m/(sm3/d)^3
    a4 = 4.13e-12;      % in m/(sm3/d)^4  
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % coefficients for a 4th order polynomial eq. (quartic equation)  %%
    % a4p*x^4 + a3p*x^3 + a2p*x^2 + a1p*x + aop, with x=fref/f        %%    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    e = a0.*nStages;
    d = a1.*(abs(qf)./(meter^3/day)).*nStages;
    c = a2.*(abs(qf)./(meter^3/day)).*nStages - dhf;
    b = a3.*(abs(qf)./(meter^3/day)).*nStages;
    a = a4.*(abs(qf)./(meter^3/day)).*nStages;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% p,q,S, and Q are terms used to solve analitically the 4th order eq %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Parameters to perform an analysis of the solutions  %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    deltaAnalysis = 256.*(a.^3).*e.^3 - 192.*(a.^2).*b.*d.*e.^2 - 128.*(a.^2).*(c.^2).*(e.^2) + 144.*(a.^2).*c.*(d.^2).*e -27.*(a.^2).*(d.^4)...
            + 144.*a.*(b.^2).*c.*(e.^2) -6.*a.*(b.^2).*(d.^2).*e - 80.*a.*b.*(c.^2).*d.*e + 18.*a.*b.*c.*(d.^3) +  16.*a.*(c.^4).*e ...
            - 4.*a.*(c.^3).*d.^2 - 27.*(b.^4).*e.^2 + 18.*(b.^3).*c.*d.*e - 4.*(b.^3).*(d.^3) - 4.*(b.^2).*(c.^3).*e + (b.^2).*(c.^2).*(d.^2);
        
    PAnalysis = 8.*a.*c - 3.*(b.^2);
    QAnalysis = b.^3 + 8.*d.*(a.^2);
    DAnalysis = 64.*(a.^3) - 16.*(a.^2).*(c.^2) + 16.*(a.^2).*b.*d - 3.*(b.^4);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Parameters to determine the solutions for the equation %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    delta0 = c.^2 - 3.*b.*d + 12.*a.*e;
    delta1 = 2.*c.^3 - 9.*b.*c.*d + 27.*(b.^2).*e + 27.*a.*d.^2 -72.*a.*c.*e;
    
    p = (8.*a.*c -3.*b.^2)./(8.*a.^2);
    q = (b.^3 - 4.*a.*b.*c + 8.*(a.^2).*d)./(8.*a.^3);
    
    Q =  ((delta1 + (delta1.^2 -4.*delta0.^3).^(1/2))./(2)).^(1/3);    
    S = 0.5.*(-(2/3).*p + (1./(3.*a)).*(Q + delta0./Q)).^(1/2);
    
    
    x1 = -b./(4.*a) - S + 0.5.*(-4.*(S.^2) - 2.*p + q./S).^(1/2);
    x2 = -b./(4.*a) - S - 0.5.*(-4.*(S.^2) - 2.*p + q./S).^(1/2);
    x3 = -b./(4.*a) + S + 0.5.*(-4.*(S.^2) - 2.*p - q./S).^(1/2);
    x4 = -b./(4.*a) + S - 0.5.*(-4.*(S.^2) - 2.*p - q./S).^(1/2);
       
    x = x1;
    
    freq = fref./x;
end

function F = root3d(x, qf, dhf)
% x = [freq; q60, dh60];
    a0 = 19.37;       % m 
    a1 = -1.23e-02;     % in m/sm3/d
    a2 = 2.24e-05;      % in m/(sm3/d)^2      
    a3 = -1.86e-08;     % in m/(sm3/d)^3
    a4 = 4.13e-12;      % in m/(sm3/d)^4
    
    qf_metric = qf./(meter^3/day);
    

    F(1) = x(1).*x(2)./(meter^3/day) - qf_metric.*60;
    F(2) = x(3) - a4.*x(2)./(meter^3/day).^4 + a3.*x(2)./(meter^3/day).^3 + a2.*x(2)./(meter^3/day).^2 + a1.*x(2)./(meter^3/day) + a0;
    F(3) = x(3).*(x(1)./60).^2 - dhf;
end