function output = from_gas_grav_to_molar_mass(input)
% output = from_gas_grav_to_molar_mass(input)
%
% Performs conversion from specific gas gravity (density relative to density of air at standard
% conditions = 1.23 kg/m^3 = 76.3 * 10^-3 lbm/ft^3) to molar mass expressed in kg/mol.
% 
% JDJ, 04-10-02, last revised 01-01-14.

output =  input * 28.97e-3;