function output = from_m2_per_d_Pa_to_bpd_psi_ft(input)
% output = from_m2_per_d_Pa_to_bpd_psi_ft(input)
%
% Performs unit conversion from m^3/(d*Pa*m) = m^2/(d*Pa) to bbl/(d*psi*ft).
%
% JDJ, 04-12-01.

output =  input / 7.565341e-5;