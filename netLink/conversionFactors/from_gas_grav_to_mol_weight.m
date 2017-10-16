function output = from_gas_grav_to_mol_weight(input)
% output = from_gas_grav_to_mol_weight(input)
%
% Performs conversion from specific gas gravity (density relative to density of air at standard
% conditions = 1.23 kg/m^3 = 76.3 * 10^-3 lbm/ft^3) to molecular weight expressed in g/mol
% (lbm/lbm-mole). 
% 
% JDJ, 01-01-14.

output =  input * 28.97;