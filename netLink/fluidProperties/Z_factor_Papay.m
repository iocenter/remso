function Z = Z_factor_Papay(p_pr,T_pr)
% Z = Z_factor_Papay(p_pr,T_pr)
%
% Calculates the Z-factor for a given pseudo reduced pressure and pseudo reduced temperature.
% Use is made of the correlation of Papay, quoted in Takacs (1976), to approximate the 
% Standing & Katz (1942) chart.
%
% Note: this is a very crude approximation to the results of Standing and Katz. A much better
% approximation is obtained using the correlation of Dranchuk & Abu-Kasem (1975), as programmed
% in Z_factor_DAK.m
%
% Z = Z factor, -
% p_pr = pseudo-reduced pressure, -
% T_pr = pseudo-reduced temperature, -
%
% JDJ, 18-09-09, last revision revion 18-09-09

Z = 1 - 3.52.*p_pr./(T_pr.*10^0.9813) + 0.274.*p_pr.^2./(T_pr.*10.^0.8157); 
