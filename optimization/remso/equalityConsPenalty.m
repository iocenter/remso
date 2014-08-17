function [rho,errorSumB,dualApproxB] = equalityConsPenalty(firstOptDualApprox,errorSum,rhoB,rhoHat,errorSumB,dualApproxB)

%L. T. Biegler, T. F. Coleman, A. R. Conn, and F. N. Santosa, Eds.,
%Large-Scale Optimization with Applications.
%Part II: Optimal Design and Control,
%vol. 93. New York, NY: Springer New York, 1997.
 


% Following suggestion by Leineweber to reduce the weights to improve
% convergence.

if errorSum < eps  %% if it is zero
    dualApprox = inf;
else
    dualApprox = abs(firstOptDualApprox)/errorSum;
end

if rhoB*errorSum >= abs(firstOptDualApprox) + 2*rhoHat*errorSum
    rho = rhoB;
elseif isinf(dualApprox)
    rho = rhoHat;
else
    rho  = dualApprox + 3*rhoHat;
end


if ~isempty(errorSumB) && isfinite(dualApprox)

    if  (dualApprox < dualApproxB) && (errorSum < errorSumB)
        rho = dualApprox + 3*rhoHat;       
    end 
    
end


errorSumB = errorSum;
dualApproxB = dualApprox;


end
