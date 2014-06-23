function [ stateG,Jac] = stdVols2stateMRST(vW,vO,vG,f,system,varargin)

opt = struct('stateG',[],'maxIter',20,'tol',1e-7);

opt = merge_options(opt, varargin{:});

if isempty(opt.stateG)  %TODO: do something smart here, use correlations
    nc = numel(vO);
    
    
    p0    =  300;
    
    p0  = convertFrom(p0, barsa);
    p0  = repmat(p0, [nc, 1]);
    s0  = repmat([ 0.2, 0.8, 0.0 ], [nc, 1]);
    rs0 = repmat( 220 , [nc, 1]);
    rv0 = 0; % dry gas
    
    stateG = struct('s', s0, 'rs', rs0, 'rv', rv0, 'pressure', p0);
    
    
else
    stateG = opt.stateG;
end

converged = false;
for  k = 1:opt.maxIter
    
    [vWG,vOG,vGG]  = blackOil2stdVols(stateG,f,system,'partials',true);
    
    eqs{1} = vWG - vW;
    eqs{2} = vOG - vO;
    eqs{3} = vGG - vG;
    
    
    [dx, linsolver_diverged,Jac] = SolveEqsADI(eqs,[]);
    if linsolver_diverged
        warning('State conversion linear solver diverged')
    end
    
    if  (norm(eqs{1}.val,'inf') < opt.tol) && (norm(eqs{2}.val,'inf') < opt.tol) && (norm(eqs{3}.val,'inf') < opt.tol)
        converged = true;
        break;
    end
    
    [stateG] = updateStateVO(stateG, dx, f, system,'updateWellSol',false);
    

end
if ~converged
    warning('State conversion did not converge')
end











end

