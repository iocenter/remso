function output = from_m2_per_s_Pa_to_bpd_psi_ft(input)
% output = from_m2_per_s_Pa_to_bpd_psi_ft(input)
%
% Performs unit conversion from m^3/(s*Pa*m) = m^2/(s*Pa) to bbl/(d*psi*ft).
%
% JDJ, 04-12-01.

output =  input /  8.756182e-10;