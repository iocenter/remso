function [errorMax, l2ErrorJo, l2ErrorJw, l2ErrorJg, l2ErrorJp] = testRunNetworkADI(E, qoE, qwE, qgE, pout, varargin)
%  TEST runNetworkADI against finite differences. 

    opt = struct('qopert',1e-06, 'qwpert', 1e-06 , 'qgpert', 1e-06, 'poutPert', 1e-06);
    opt = merge_options(opt, varargin{:});    
   
    
   [qoE, qwE, qgE, pout] = initVariablesADI(qoE, qwE, qgE, pout);    
    
    
    dp = dpBeggsBrill(E, qoE, qwE, qgE, pout);  

%     dp = simpleDp(E,qoE, qwE, qgE, pout);
    
    jo = cell2mat(dp.jac(1)); % oil rate jacobian
    jw = cell2mat(dp.jac(2)); % water rate jacobian
    jg = cell2mat(dp.jac(3)); % gas rate jacobian   
    jp = cell2mat(dp.jac(4)); % TODO: evaluate the part of the jacobian related to average pressure

    
    fo = @(qoE) dpBeggsBrill(E, qoE, qwE, qgE, pout);
    fw = @(qwE) dpBeggsBrill(E, qoE, qwE, qgE, pout); 
    fg = @(qgE) dpBeggsBrill(E, qoE, qwE, qgE, pout); 
    fp = @(pout) dpBeggsBrill(E, qoE, qwE, qgE, pout); 
    
    [dfdqo] = calcPertGrad(fo,qoE.val,opt.qopert);    
    [dfdqw] = calcPertGrad(fw,qwE.val,opt.qwpert);    
    [dfdqg] = calcPertGrad(fg,qgE.val,opt.qgpert);
    [dfdp] = calcPertGrad(fp,pout.val,opt.poutPert);

%     l2ErrorJo = abs((dfdqo-jo)/dfdqo);
%     l2ErrorJw = abs((dfdqw-jw)/dfdqw);    
%     l2ErrorJg = abs((dfdqg-jg)/dfdqg);
    
    l2ErrorJo = abs(dfdqo-jo)
    jo
    l2ErrorJw = abs(dfdqw-jw)
    jw
    l2ErrorJg = abs(dfdqg-jg)
    jg
    l2ErrorJp = abs(dfdp-jp)
    jp
    
    errorMax = max(max(l2ErrorJo,l2ErrorJw), max(l2ErrorJg, l2ErrorJp));

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

