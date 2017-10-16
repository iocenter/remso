function output = from_Pa_per_m_to_kg_per_m3(input)
% output = from_Pa_per_m_to_kg_per_m3(input)
%
% Performs conversion from gradient in Pa/m
% to density in kg/m^3
%
% JDJ, 29-11-01.

output =  input / 9.80665;