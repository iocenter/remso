function output = from_bpd_psi_ft_to_m2_per_d_Pa(input)
% function output = from_bpd_psi_ft_to_m2_per_d_Pa(input)
%
% Performs unit conversion from bbl/(d*psi*ft) to m^3/(d*Pa*m) = m^2/(d*Pa).
%
% JDJ, 19-11-01, last revised 10-05-13

output = input * 7.565341e-5;