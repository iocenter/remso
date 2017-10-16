function output = from_lbm_per_ft3_to_Pa_per_m(input)
% output = from_lbm_per_ft3_to_Pa_per_m(input)
%
% Performs unit conversion from density in lbm/ft^3 to
% gradient in Pa/m.
%
% JDJ, 29-11-01.

output =  input * 1.601846e1 * 9.80665 ;