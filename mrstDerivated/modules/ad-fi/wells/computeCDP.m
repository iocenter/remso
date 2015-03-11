function [ wellSol ] = computeCDP( W, wellSol, bhp, q_s, p, b, r, m,rMax, rho_s, model,opt,varargin )

opt.maxIts = 25;
opt.tol = 0.01*Pascal;

opt = merge_options(opt, varargin{:});


%instead of fixing cdp, lets calculate it
if strcmp(opt.cdpCalc,'exact')
    % clean all jacobians of the primary variables
    toDouble = @(x)cellfun(@double, x, 'UniformOutput', false);
    wv.q_s = toDouble(q_s);
    wv.bhp = double(bhp);
    wv.b = toDouble(b);
    wv.m = toDouble(m);
    wv.p = double(p);
    wv.rMax = toDouble(rMax);
    wv.r = toDouble(r);
    
    
    nConn       = cellfun(@numel, {W.cells})'; % # connections of each well
    perf2well   = rldecode((1:numel(W))', nConn);
    ixW = arrayfun(@(wnr)perf2well==wnr,1:numel(wellSol),'UniformOutput',false);


    % solve for cdp first (and consequently for cqs).
    % Solve f(cdp,p) - cdp = 0  for a fixed p (p stands for primary variables, the outer loop variables)
    wellSol = arrayfun(@(wsi)subsasgn(wsi,struct('type','.','subs','cdp'),double(wsi.cdp)),wellSol);
    nW = numel(wellSol);
    its = 1 ;
    converged = false;
    while ~converged && (opt.maxIts >= its)
        
        % assume cdp is independent
        [wellSol.cdp] = initVariablesADI(wellSol.cdp);
        cdp = {wellSol.cdp};
        
        
        [~, cq_s, ~, ~, ~, ~] = computeWellContributions(...
            W, wellSol, wv.bhp, wv.q_s, wv.p, wv.b, wv.r, wv.m, model, ...
            opt.allowWellSignChange, opt.allowCrossFlow);
        
        wellSol = arrayfun(@(wsi,wn)subsasgn(wsi,struct('type',{'.'},'subs',{'cqs'}),...    wellSol(i).cqs = cq_s(ixW(i))
                                   cellfun(@(ci,ix)ci(ixW{wn}),cq_s,'UniformOutput',false)),...
                            wellSol,(1:numel(wellSol))');
        
        wellSol = updateConnDP(W, wellSol, wv.b, wv.rMax, rho_s, model);
        
        % This is  f(cdp) - cdp  (which we want to take to be == 0)
        cdpDiff = cellfun(@(w,c)w-c,{wellSol.cdp},cdp,'UniformOutput',false);
        
        %Compute the Newton correction
        %cdpDiff is block diagonal. Exploit it!
        deltaCdp = cellfun(@(c,k)-c.jac{k}\c.val,cdpDiff,num2cell(1:nW),'UniformOutput',false);
        
        % Add correction
        cdp = cellfun(@(c,d)double(c)+d,cdp,deltaCdp,'UniformOutput',false);
        [wellSol.cdp] = cdp{:};
        
        % check convergence
        residual = norm(cell2mat(deltaCdp));
        converged = residual < opt.tol;
        
        its = its +1;
    end
    
    if converged && isa(bhp,'ADI')  %% correct the jacobians.  Implicit function solved!
        % Compute \frac{d f}{d p} i.e.  The total derivative of
        % f w.r.t to the primary variables
        
        % Compute \frac{\partial f}{\partial p}  , assuming cdp independent
        [~, cq_s, ~, ~, ~, ~] = computeWellContributions(...
            W, wellSol, bhp, q_s, p, b, r, m, model, ...
            opt.allowWellSignChange, opt.allowCrossFlow);
        
        wellSol = arrayfun(@(wsi,wn)subsasgn(wsi,struct('type',{'.'},'subs',{'cqs'}),...    wellSol(i).cqs = cq_s(ixW(i))
                                   cellfun(@(ci,ix)ci(ixW{wn}),cq_s,'UniformOutput',false)),...
                            wellSol,(1:numel(wellSol))');
        
        wellSol = updateConnDP(W, wellSol, b, rMax, rho_s, model);
        
        
        
        % Compute df/dp = - (\frac{\partial f}{\partial d} - I )^{-1}* \frac{\partial f}{\partial p}
        
        invDeltaCdp = cellfun(@(c,k)-inv(c.jac{k}),cdpDiff,num2cell(1:nW),'UniformOutput',false);
        wellSol =   arrayfun(@(ws,k)subsasgn(ws,struct('type',{'.','.'},'subs',{'cdp','jac'}),...  wellSol(w).cdp.jac(:) =
            cellfun(@(j) invDeltaCdp{k}*j,ws.cdp.jac,'UniformOutput',false)),wellSol,(1:numel(wellSol))');  % invDeltaCdp(w) * wellSol(w).cdp.jac(:)
        % now the jacobian of cdp w.r.t the primary variables is finally correct!,
        % compute cqs with correct cdp and jacobians
        
        
    elseif ~converged
        % Things we
        warning(['Well equations did not converge. norm(residual) = ' num2str(residual) ' Pascal']);
        
        
        % Update well pressure according to MRST default
        if opt.iteration ==1
            [wellSol, ~, ~] = updateConnDP(W, sol, b, rMax, rhos, model);
        end
        
    end
elseif strcmp(opt.cdpCalc,'none')
    % do nothing
elseif strcmp(opt.cdpCalc,'first')
    if opt.iteration ==1    
        [wellSol, ~, ~] = updateConnDP(W, sol, b, rMax, rhos, model);
    end
else
	error(['Could not recognize cdpCalc method. cdpCalc == ' opt.cdpCalc]);
end




end

