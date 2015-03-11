function [ wellSol,wellmodel ] = computeCDPOOP( wellmodel,currentFluxes,bhp,wellSol,model )


%instead of fixing cdp, lets calculate it
if strcmp(wellmodel.cdpCalc,'exact')
    % clean all jacobians of the primary variables
    toDouble = @(x)cellfun(@double, x, 'UniformOutput', false);
    currentFluxesVal = toDouble(currentFluxes);
    bhpVal = double(bhp);
    wellmodelVal = wellmodel;
    wellmodelVal.bfactors = toDouble(wellmodel.bfactors);
    wellmodelVal.mobilities = toDouble(wellmodel.mobilities);
    wellmodelVal.pressure = toDouble(wellmodel.pressure);
    wellmodelVal.referencePressure = double(wellmodel.referencePressure);
    wellmodelVal.maxComponents = toDouble(wellmodel.maxComponents);
    wellmodelVal.components = toDouble(wellmodel.components);
    wellmodelVal.saturations = toDouble(wellmodel.saturations);
    
    
    % solve for cdp first (and consequently for cqs).
    % Solve f(cdp,p) - cdp = 0  for a fixed p (p stands for primary variables, the outer loop variables)
    wellSol = arrayfun(@(wsi)subsasgn(wsi,struct('type','.','subs','cdp'),double(wsi.cdp)),wellSol);
    nW = numel(wellSol);
    its = 1 ;
    converged = false;
    while ~converged && (wellmodel.maxIts >= its)
        
        % assume cdp is independent
        [wellSol.cdp] = initVariablesADI(wellSol.cdp);
        cdp = {wellSol.cdp};
        
        
        [~, ~, ~, ~, wellSol] =...
            wellmodelVal.assembleEquations(wellSol, currentFluxesVal, bhpVal, model);
        [wellSol, ~, ~] = wellmodelVal.updatePressure(wellSol, currentFluxesVal, bhpVal, model);
        
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
        converged = residual < wellmodel.tol;
        
        its = its +1;
    end
    
    if converged && isa(bhp,'ADI')  %% correct the jacobians.  Implicit function solved!
        % Compute \frac{d f}{d p} i.e.  The total derivative of
        % f w.r.t to the primary variables
        
        % Compute \frac{\partial f}{\partial p}  , assuming cdp independent
        [~, ~, ~, ~, wellSol] =...
            wellmodelVal.assembleEquations(wellSol, currentFluxes, bhp, model);
        [wellSol, ~, ~] = wellmodel.updatePressure(wellSol, currentFluxes, bhp, model);
        
        
        % Compute df/dp = - (\frac{\partial f}{\partial d} - I )^{-1}* \frac{\partial f}{\partial p}
        
        invDeltaCdp = cellfun(@(c,k)-inv(c.jac{k}),cdpDiff,num2cell(1:nW),'UniformOutput',false);
        wellSol =   arrayfun(@(ws,k)subsasgn(ws,struct('type',{'.','.'},'subs',{'cdp','jac'}),...  wellSol(w).cdp.jac(:) =
            cellfun(@(j) invDeltaCdp{k}*j,ws.cdp.jac,'UniformOutput',false)),wellSol,1:numel(wellSol));  % invDeltaCdp(w) * wellSol(w).cdp.jac(:)
        % now the jacobian of cdp w.r.t the primary variables is finally correct!,
        % compute cqs with correct cdp and jacobians
        
        
    elseif ~converged
        % Things we
        warning(['Well equations did not converge. norm(residual) = ' num2str(residual) ' Pascal']);
        
        
        % Update well pressure according to MRST default
        [wellSol, ~, ~] = wellmodel.updatePressure(wellmodel, wellSol, currentFluxes, bhp, model);
    end
elseif strcmp(wellmodel.cdpCalc,'none')
    % do nothing
elseif strcmp(wellmodel.cdpCalc,'fixed')
        [wellSol, ~, ~] = wellmodel.updatePressure(wellmodel, wellSol, currentFluxes, bhp, model);
else
	error(['Could not recognize cdpCalc method. cdpCalc == ' wellmodel.cdpCalc]);
end




end

