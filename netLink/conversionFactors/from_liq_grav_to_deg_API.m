function output = from_liq_grav_to_deg_API(input)
% output = from_liq_grav_to_deg_API(input)
%
% Performs unit conversion from liquid gravity (wrt. water = 1000 kg/m^3) to degrees API.
%
% JDJ, 03-02-03.

output =  141.5e3/(1000*input) - 131.5;