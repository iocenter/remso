function output = from_bpd_psi_ft_to_m2_per_s_Pa(input)
% function output = from_bpd_psi_ft_to_m2_per_s_Pa(input)
%
% Performs unit conversion from bbl/(d*psi*ft) to m^3/(s*Pa*m) = m^2/(s*Pa).
%
% JDJ, 19-11-01, last revised 10-05-13

output = input * 8.756182e-10;