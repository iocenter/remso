function p_b = pres_bub_Standing(R_sb,rho_g_sc,rho_o_sc,T)
% function p_b = pres_bub_Standing(R_sb,rho_g_sc,rho_o_sc,T)
%
% Computes the bubble point pressure with a Standing correlation converted to SI units.
%
% R_sb = gas-oil ratio at bubble point pressure, m^3/m^3
% p_b = bubble point pressure, Pa
% rho_g_sc = gas density at standard conditions, kg/m3
% rho_o_sc = oil density at standard conditions, kg/m3
% T = temperature, deg. C
%
% JDJ, 27-02-01, last revised 10-05-13
%{
%Changes by Thiago and Codas
%Make the function compatible with ADI objects
%Vectorization
%}
% check for presence of gas:
p_b = R_sb; % initialize variable
c1 = (rho_g_sc == 0 | (double(R_sb) == 0));
p_b(c1) =  1.e5; % atmospheric pressure
    
if any(~c1)
    help01 = (10.^(0.00164*T))/(10.^(1768./rho_o_sc));
    help01 = help01.*ones(size(double(R_sb)));
    p_b(~c1) = 125e3 * ((716.*R_sb(~c1)./rho_g_sc).^0.83 .* help01(~c1) - 1.4);
end

% Reality check:
p_b(p_b < 1.e5) = 1.e5; % atmospheric pressure

end
