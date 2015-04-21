function [ error ] = testSingleShooting(u,ss,obj,varargin)


opt = struct('pert',1e-6,'debug',false);
opt = merge_options(opt, varargin{:});

uDims = cellfun(@numel,u);

[f,gradU,converged,simVarsOut,xs,vs] = simulateSystemSS(u,ss,obj,'gradients',true);

fu = @(uu) simulateSystemSS(mat2cell(uu,uDims,1),ss,obj,'gradients',false);
dfdu = calcPertGrad(fu,cell2mat(u),opt.pert);


uRightSeeds = speye(sum(uDims));
uRightSeeds = mat2cell(uRightSeeds,uDims,sum(uDims));


[f,gradUF] = simulateSystemSS(u,ss,obj,'gradients',true,'guessV',vs,'guessX',xs,'simVars',simVarsOut,'uRightSeeds',uRightSeeds);



error = max([norm(cell2mat(gradU) - dfdu),norm(cell2mat(gradU)-gradUF)]);

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