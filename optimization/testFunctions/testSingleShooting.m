function [ error ] = testSingleShooting(u,ss,obj,varargin)


opt = struct('pert',1e-6,'debug',false);
opt = merge_options(opt, varargin{:});

nu = numel(u{1});

[f,gradU,converged,simVarsOut] = simulateSystemSS(u,ss,obj,'gradients',true);

fu = @(uu) simulateSystemSS(toStructuredCells(uu,nu),ss,obj,'gradients',false);
dfdu = calcPertGrad(fu,cell2mat(u),opt.pert);


error = norm(cell2mat(gradU) - dfdu);

end


function [dfdx] = calcPertGrad(f,xb,pert)

nx = numel(xb);
fb = f(xb);

dfdx = zeros(numel(fb),nx);

parfor j = 1:nx
    xp = xb;
    xp(j) = xp(j) + pert;
    
    fp = f(xp);
    
    dfdx(:,j) = (fp-fb)/pert;
    
end

end