function [freq, q60, dh60 ] = pump_eq_system(qf, dhf)
% Solves a system of equations to obtain frequency, flow and dh at 60Hz
% x = [freq; q60, dh60];

    fun = @root3d;    
    options = optimoptions('fsolve','Display','iter');
    
    freq = zeros(numel(qf),1);
    q60 = zeros(numel(qf),1);
    dh60 = zeros(numel(qf),1);
    x = zeros(numel(qf),3);
    for i=1:numel(qf)        
        x0 = [50; qf(i); dhf(i)];
        
        x(i,:) = fsolve(@(x) fun(x, qf(i), dhf(i)), x0);
    end
    
    freq = x(:,1);
    q60 = x(:,2);
    dh60 = x(:,3);
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