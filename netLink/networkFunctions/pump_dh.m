function [ dhf ] = pump_dh(dpf, rho)
% Calculates the head difference for the pump at a given pressure loss.
% It can be used to compute upper and lower bounds for the pressure loss
% in the equipment. The coefficients of the curve were found in the paper
% 'Exploring the potential of model-based optimization in oil production
% gathering networks with esp-produced, high water cut wells'.
%    
    dhf = abs(dpf)./(rho.*norm(gravity));
end



