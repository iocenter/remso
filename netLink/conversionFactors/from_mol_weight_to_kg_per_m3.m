function output = from_mol_weight_to_kg_per_m3(input)
% output = from_mol_weight_to_kg_per_m3(input)
%
% Performs conversion from molecular weight expressed in g/mol (lbm/lbm-mole) to kg/m^3.
%  
% JDJ, 01-01-14

output =  input / 23.55;