function [tau] = boundViolationWeights(mu,tauB,withAlgs)
%  Calculation of the penalty parameter (tau) for the inequality
%  constraints according to:
%  
%  M. Powell, “A fast algorithm for nonlinearly constrained optimization
%              calculations,” Numer. Anal., vol. 630, pp. 144–157, 1978.
%  
%  Equation (4.6) of the paper
%
%  tau_it = max(mu,(mu+tau_{it-1})/2);
%
if withAlgs
    tauN = {...
        cellfun(@(x1,x2)abs(x1-x2),mu.ub.x,mu.lb.x,'UniformOutput',false);...
        cellfun(@(x1,x2)abs(x1-x2),mu.ub.v,mu.lb.v,'UniformOutput',false)};
else
    tauN = ...
        {cellfun(@(x1,x2)abs(x1-x2),mu.ub.x,mu.lb.x,'UniformOutput',false)};
end

if isempty(tauB)
    tau = tauN;
else
    tau = cellfun(@(tauV,tauVB) cellfun(@(x1,x2)max(x1,(x1+x2)/2),tauV,tauVB,'UniformOutput',false),tauN,tauB,'UniformOutput',false);
end

end