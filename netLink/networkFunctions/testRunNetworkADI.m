function [errorMax, joRelError, jwRelError, jgRelError, jpReError] = testRunNetworkADI(E, qoE, qwE, qgE, pout, varargin)
%  TEST runNetworkADI against finite differences. 

    opt = struct('qopert',1e-06, 'qwpert', 1e-06 , 'qgpert', 1e-06, 'poutPert', 1e-06,'dpFunction', @dpBeggsBrill);
    opt = merge_options(opt, varargin{:});    
   
    
   [qoE, qwE, qgE, pout] = initVariablesADI(qoE, qwE, qgE, pout);    
    
    
    dp = opt.dpFunction(E, qoE, qwE, qgE, pout);  

%     dp = simpleDp(E,qoE, qwE, qgE, pout);
    
    jo = cell2mat(dp.jac(1)); % oil rate jacobian
    jw = cell2mat(dp.jac(2)); % water rate jacobian
    jg = cell2mat(dp.jac(3)); % gas rate jacobian   
    jp = cell2mat(dp.jac(4)); % average pressure jacobian

    
    fo = @(qoE)  opt.dpFunction(E, double(qoE), double(qwE), double(qgE), double(pout));
    fw = @(qwE)  opt.dpFunction(E, double(qoE), double(qwE), double(qgE), double(pout)); 
    fg = @(qgE)  opt.dpFunction(E, double(qoE), double(qwE), double(qgE), double(pout)); 
    fp = @(pout) opt.dpFunction(E, double(qoE), double(qwE), double(qgE), double(pout)); 
    
    [dfdqo] = calcPertGrad(fo,qoE.val,opt.qopert);    
    [dfdqw] = calcPertGrad(fw,qwE.val,opt.qwpert);    
    [dfdqg] = calcPertGrad(fg,qgE.val,opt.qgpert);
    [dfdp] = calcPertGrad(fp,pout.val,opt.poutPert);

%     l2ErrorJo = abs((dfdqo-jo)/dfdqo);
%     l2ErrorJw = abs((dfdqw-jw)/dfdqw);    
%     l2ErrorJg = abs((dfdqg-jg)/dfdqg);
    
    l2ErrorJo = abs(dfdqo-jo);
    jo, dfdqo
    l2ErrorJw = abs(dfdqw-jw);
    jw, dfdqw
    l2ErrorJg = abs(dfdqg-jg);
    jg, dfdqg
    l2ErrorJp = abs(dfdp-jp);
    jp, dfdp

    if jo ~= 0        
        joRelError = norm(l2ErrorJo/jo);
    else 
        joRelError = l2ErrorJo;
    end
    
    if jw ~=0
        jwRelError = norm(l2ErrorJw/jw);
    else 
        jwRelError = l2ErrorJw;
    end
    
    if jg ~=0
        jgRelError = norm(l2ErrorJg/jg);
    else
        jgRelError = l2ErrorJg;
    end
    
    if jp ~= 0
        jpReError = norm(l2ErrorJp/jp);        
    else        
        jpReError = l2ErrorJp;
    end    
   
    errorMax = max(max(joRelError,jwRelError), max(jgRelError, jpReError));

end


function [dfdx] = calcPertGrad(f,xb,pert)

    nx = numel(xb);
    fb = f(xb);

    dfdx = zeros(numel(fb),nx);

    for j = 1:nx
        xp = xb;
        xp(j) = xp(j) + pert;

        fp = f(xp);

        dfdx(:,j) = (fp-fb)/pert;

    end

end

