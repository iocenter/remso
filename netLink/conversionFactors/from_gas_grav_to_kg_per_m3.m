function output = from_gas_grav_to_kg_per_m3(input)
% output = from_gas_grav_to_kg_per_m3(input)
%
% Performs conversion from specific gas gravity (density relative
% to density of air at standard conditions = 1.23 kg/m^3 =
% 76.3 * 10^-3 lbm/ft^3) to kg/m^3.
% 
% JDJ, 30-10-01.

output =  input * 1.23;