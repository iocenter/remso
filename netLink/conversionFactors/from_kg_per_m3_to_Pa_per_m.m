function output = from_kg_per_m3_to_Pa_per_m(input)
% output = from_kg_per_m3_to_Pa_per_m(input)
%
% Performs conversion from density in kg/m^3 to
% gradient in Pa/m.
%
% JDJ, 29-11-01.

output =  input * 9.80665;