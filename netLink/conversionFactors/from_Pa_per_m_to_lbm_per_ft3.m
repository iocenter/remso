function output = from_Pa_per_m_to_lbm_per_ft3(input)
% output = from_Pa_per_m_to_lbm_per_ft3(input)
%
% Performs unit conversion from gradient in Pa/m to
% density in lbm/ft^3.
%
% JDJ, 29-11-01.

output =  input / (9.80665 * 1.601846e1);