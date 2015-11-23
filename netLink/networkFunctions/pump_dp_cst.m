function [ dp ] = pump_dp_cst(qp, rho)
% Calculates the pressure loss in the pump for a given flow.
% It can be used to compute upper and lower bounds for the pressure loss
% in the equipment. The coefficients of the curve were found in the paper
% 'Exploring the potential of model-based optimization in oil production
% gathering networks with esp-produced, high water cut wells'.
%    

    dh = pump_dh(qp); % in m    
    
    dp = rho.*norm(gravity).*dh;
end

