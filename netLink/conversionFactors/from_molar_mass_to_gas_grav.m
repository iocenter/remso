function output = from_molar_mass_to_gas_grav(input)
% output = from_molar_mass_to_gas_grav(input)
%
% Performs conversion from molar mass expressed in kg/mol to specific gas gravity (density
% relative to density of air at standard conditions = 1.23 kg/m^3 = 76.3 * 10^-3 lbm/ft^3).
% 
% JDJ, 04-10-02, last revised 01-01-14.

output =  input / 28.97e-3;