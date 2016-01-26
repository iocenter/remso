function sigma = interfacial_tensions()
% function sigma = interfacial_tensions()
%
% Input function for interfacial tensions
%
% JDJ, 17-03-03, last revised 10-05-13 

sigma_go = 0.008; % gas-oil interfacial tension, N/m
sigma_gw = 0.04; % gas-water interfacial tension, N/m
sigma = [sigma_go,sigma_gw];
