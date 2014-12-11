function [rho,errorSumB,dualApproxB] = equalityConsPenalty(gbarR,errorSum,crossProduct,rhoB,rhoHat,errorSumB,dualApproxB)
%
%
% dualApprox = (abs((g + nu)' Y p_y) + abs(min(w'*pz/(1-xi),0)))  / norm(c,1)
%


dualApprox = (abs(gbarR)+abs(min(crossProduct,0)))/errorSum;

if ~isfinite(dualApprox)
    rho = max(rhoHat,rhoB);
elseif rhoB*errorSum >= dualApprox + 2*rhoHat*errorSum
    rho = rhoB;
else
    rho  = dualApprox + 3*rhoHat;
end



% Following suggestion by Leineweber to reduce the weights to improve
% convergence.


%
% D. B. Leineweber,
% "Efficient Reduced SQP Methods for the Optimization of Chemical Processes Described by Large Sparse DAE Models," 
% University of Heidelberg, 
% 1998.
%
% Page 129, suggestion to decrease weights
%

if ~isempty(errorSumB) && isfinite(dualApprox)

    if  (dualApprox < dualApproxB) && (errorSum < errorSumB)
        rho = dualApprox + 3*rhoHat;       
    end 
    
end


errorSumB = errorSum;
dualApproxB = dualApprox;


end
