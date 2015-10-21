function [errorMax, l2ErrorJo, l2ErrorJw, l2ErrorJg] = testRunNetworkADI(E, qoE, qwE, qgE, pout, varargin)
%  TEST runNetworkADI against finite differences. 

    opt = struct('qopert',1e-06, 'qwpert', 1e-06 , 'qgpert', 1e-06);
    opt = merge_options(opt, varargin{:});    
   
    
   [qoE, qwE, qgE, pout] = initVariablesADI(qoE, qwE, qgE, pout);
 
    
    
    
%     dp = dpBeggsBrill(E, qoE, qwE, qgE, pout);  

    dp = simpleDp(E,qoE, qwE, qgE, pout);
    
    jo = cell2mat(dp.jac(1:length(qoE.val)));
    jw = cell2mat(dp.jac(length(qoE.val)+1:length(qoE.val)+length(qwE.val)));
    jg = cell2mat(dp.jac(length(qoE.val)+length(qwE.val)+1:length(qoE.val)+length(qwE.val)+length(qgE.val)));
    
    fo = @(qoE) simpleDp(E, qoE, qwE, qgE, pout);
    fw = @(qwE) simpleDp(E, qoE, qwE, qgE, pout); 
    fg = @(qgE) simpleDp(E, qoE, qwE, qgE, pout); 
    
    
    [dfdqo] = calcPertGrad(fo,qoE.val,opt.qopert);    
    [dfdqw] = calcPertGrad(fw,qwE.val,opt.qwpert);    
    [dfdqg] = calcPertGrad(fg,qgE.val,opt.qgpert);

    
    l2ErrorJo = norm(abs((dfdqo-jo)));
    
    l2ErrorJw = norm(abs((dfdqw-jw)));
    
    l2ErrorJg = norm(abs((dfdqg-jg)));
    
    errorMax = max(max(l2ErrorJo,l2ErrorJw), l2ErrorJg);

end


function [dfdx] = calcPertGrad(f,xb,pert)

    nx = numel(xb);
    fb = f(xb);

    dfdx = zeros(numel(fb.val),nx);

    for j = 1:nx
        xp = xb;
        xp(j) = xp(j) + pert;

        fp = f(xp);

        dfdx(:,j) = (fp.val-fb.val)/pert;

    end

end

