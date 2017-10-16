function output = from_mol_weight_to_gas_grav(input)
% output = from_mol_weight_to_gas_grav(input)
%
% Performs conversion from molecular weight expressed in g/mol (lbm/lbm-mole) to specific gas
% gravity (density relative to density of air at standard conditions = 1.23 kg/m^3 = 76.3 *
% 10^-3 lbm/ft^3). 
% 
% JDJ, 01-01-14.

output =  input / 28.97;