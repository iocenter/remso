function output = from_deg_API_to_liq_grav(input)
% function output = from_deg_API_to_liq_grav(input)
%
% Performs unit conversion from degrees API to liquid gravity (wrt. water = 1000 kg/m^3).
%
% JDJ, 03-02-03, last revised 10-05-13

output = (141.5e3/(131.5+input)) / 1000;