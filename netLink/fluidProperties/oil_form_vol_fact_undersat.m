function B_o = oil_form_vol_fact_undersat(B_ob,c_o,p,p_b)
% function B_o = oil_form_vol_fact_undersat(B_ob,c_o,p,p_b)
%
% Computes the oil formation volume factor for undersaturated oil. Valid for SI units and field
% units. 
%
% B_o = oil formation volume factor, m^3/m^3, (bbl/bbl)
% B_ob = oil formation volume factor at bubble point pressure,
%        m^3/m^3, (bbl/bbl)
% c_o = compressibility, 1/Pa, (1/psi)
% p = presssure, Pa, (psi)
% p_b = bubble point pressure, Pa, (psi)
%
% JDJ, 02-03-01, last revised 10-05-13
%
%Change by codas: Vectorization

B_o = B_ob .* exp(-c_o.*(p - p_b));