function varargout = condensing(x,u,v,ss,varargin)
% Apply the Lift-Opt trick and get the correction and predictor matrices.
%
% SYNOPSIS:
%  [xs,vs,xd,vd,ax,Ax,av,Av] = condensing(x,u,v,ss)
%  [xs,vs,xd,vd,ax,Ax,av,Av] = condensing(x,u,v,ss, 'pn', pv, ...)
%
% PARAMETERS:
%
%   x - States in the prediction horizon. Initial conditions for
%       the shooting intervals (except for the first which a known constant)
%
%   u - cellarray containing the controls for each control
%       period.
%
%   v - Algebraic states in the prediction horizon.
%
%   ss - A simulator structure, containing all the required
%        information on the model.
%
%   'pn'/pv - List of 'key'/value pairs defining optional parameters. The
%             supported options are:
%
%   simVars - Provide the simulation result to hot-start the Jacobian
%             calculation. If not given, the simulation will be
%             recalculated.  Suggestion: run symulateSystem to get it
%             before!
%
%
% RETURNS:
%
%   xs - Simulated state output.
%
%   vs - Simulated algebraic state output.
%
%   xd - Difference xs-x. (Multiple Shooting error)
%
%   vd - Difference vs-v. (Multiple Shooting error)
%
%   ax Ax - A predicition model for the states
%           \Delta x = x + ax + Ax * \Delta u
%
%   av Av - A predicition model for the algebraic states
%           \Delta v = v + av + Av * \Delta u
%
% SEE ALSO:
%
%

opt = struct('simVars',[],'withAlgs',false);
opt = merge_options(opt, varargin{:});


% Pepare inputs

totalPredictionSteps = getTotalPredictionSteps(ss);
totalControlSteps = numel(u);

nx = numel(ss.state);
uDims =  cellfun(@(z)numel(z),u);

withAlgs = opt.withAlgs;

xd = cell(totalPredictionSteps,1);
vd = cell(totalPredictionSteps,1);
xs = cell(totalPredictionSteps,1);
vs = cell(totalPredictionSteps,1);


if isempty(opt.simVars)
    simVars = cell(totalPredictionSteps,1);
elseif iscell(opt.simVars)
    simVars = opt.simVars;
else
    simVars = bringVariables(opt.simVars,ss.jobSchedule);
end

dzdd = zeros(nx,1);

ax = cell(totalPredictionSteps,1);

Ax = cell(totalPredictionSteps,totalControlSteps);

av = cell(totalPredictionSteps,1);

Av = cell(totalPredictionSteps,totalControlSteps);


converged = false(totalPredictionSteps,1);

% Gradient propagation
t0 = tic;
k0 = 0;
for k = 1:totalPredictionSteps
    [t0,k0] = printCounter(1, totalPredictionSteps, k,'Condensing',t0,k0);
    
    
    i = callArroba(ss.ci,{k});
    ui = u{i};
    
    % Prepare the RHS for the condensing technique
    
    if k >1
        if isempty(Ax{k-1,i})
            xRightSeeds = [{dzdd},Ax(k-1,1:i-1),{zeros(nx,uDims(i))}];            
        else
            xRightSeeds = [{dzdd},Ax(k-1,1:i)];
        end       
        
        seedSizes = cellfun(@(x)size(x,2),xRightSeeds);   
        xRightSeeds = cell2mat(xRightSeeds);
        
        uRightSeeds = [zeros(uDims(i),sum([1;uDims(1:i-1)])), eye(uDims(i))];

        xStart = x{k-1};
        
    elseif  k == 1
        
        xRightSeeds = zeros(nx,uDims(1));      
        uRightSeeds = eye(uDims(1));

        xStart = ss.state;
        
    else
        error('what?')
    end
    

    
    % Compute the Jacobian-vector products
    [xs{k},vs{k},Jac,convergence] = ss.stepClient{k}(xStart,ui,'gradients',true,'xRightSeeds',xRightSeeds,'uRightSeeds',uRightSeeds,'simVars',simVars{k});
    
    converged(k) = convergence.converged;
    
    xd{k} = (xs{k}-x{k});
    

    if withAlgs
        vd{k} = (vs{k}-v{k});
    end
    
    % Extract the linearized model information
    if k >1
        [dzddAx] = mat2cell(Jac.xJ,nx,seedSizes);
        [dzdd,Ax{k,1:i}] = deal(dzddAx{:});
    else % k ==1
        [Ax{1,1}] = Jac.xJ;        
    end
    
    dzdd = dzdd -xd{k};
    
    ax{k}    = -dzdd;
    
    if withAlgs
        if k >1
            nv = size(Jac.vJ,1);
            [dvddvsu] = mat2cell(Jac.vJ,nv,seedSizes);
            [dvdd,Av{k,1:i}] = deal(dvddvsu{:});
            
            av{k}    = vd{k}-dvdd;
        else %k==1
            Av{1,1} = Jac.vJ;
            av{1}    = vd{k};
        end
        
    end
    
    
end

if ~all(converged)
    steps = 1:totalPredictionSteps;
    warning(strcat('Steps failed to converge:',num2str(steps(~converged))));
end

varargout{1} = xs;
varargout{2} = vs;
varargout{3} = xd;
varargout{4} = vd;
varargout{5} = ax;
varargout{6} = Ax;
varargout{7} = av;
varargout{8} = Av;



end

