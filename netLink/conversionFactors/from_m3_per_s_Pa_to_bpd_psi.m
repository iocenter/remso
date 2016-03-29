function output = from_m3_per_s_Pa_to_bpd_psi(input)
% output = from_m3_per_s_Pa_to_bpd_psi(input)
%
% Performs unit conversion from m^3/(s*Pa) to bbl/(d*psi).
%
% JDJ, 29-11-01.

output =  input /  2.668884e-10;