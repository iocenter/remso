function output = from_deg_API_to_kg_per_m3(input)
% function output = from_deg_API_to_kg_per_m3(input)
%
% Performs unit conversion from degrees API to kg/m^3.
%
% JDJ, 29-11-01, last revised 10-05-13

output = 141.5e3/(131.5+input);