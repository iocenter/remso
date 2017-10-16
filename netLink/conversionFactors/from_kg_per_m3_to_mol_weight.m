function output = from_kg_per_m3_to_mol_weight(input)
% output = from_kg_per_m3_to_mol_weight(input)
%
% Performs conversion from kg/m^3 to molecular weight expressed in g/mol (lbm/lbm-mole).
%  
% JDJ, 01-01-14.

output =  input * 23.55;