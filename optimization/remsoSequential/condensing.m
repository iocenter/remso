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
%   uRightSeeds - return Ax*uRightSeeds and Av*uRightSeeds instead of Ax and
%                 Av 
%
%   computeCorrection - true if ax and av should be computed
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

opt = struct('simVars',[],'uRightSeeds',[],'computeCorrection',false,'computeNullSpace',true,'xd',[],'vd',[],'withAlgs',false);
opt = merge_options(opt, varargin{:});


% Pepare inputs

totalPredictionSteps = getTotalPredictionSteps(ss);
totalControlSteps = numel(u);

nx = numel(ss.state);
uDims =  cellfun(@(z)numel(z),u);

withAlgs = opt.withAlgs;

givenRangeRHS = ~isempty(opt.xd);

if givenRangeRHS
    xd = opt.xd;
    vd = opt.vd;
else
xd = cell(totalPredictionSteps,1);
vd = cell(totalPredictionSteps,1);
end
xs = cell(totalPredictionSteps,1);
vs = cell(totalPredictionSteps,1);


if isempty(opt.simVars)
    simVars = cell(totalPredictionSteps,1);
else
    simVars = opt.simVars;
end

if ~opt.computeNullSpace
    opt.uRightSeeds = cellfun(@(xx)zeros(numel(xx),0),u,'UniformOutput',false);
end

uSeedsProvided = true;
if isempty(opt.uRightSeeds)
    uSeedsProvided = false; 
    nuSeed = cellfun(@(xx)numel(xx),u);
else
    nuSeed = cellfun(@(xx)size(xx,2),opt.uRightSeeds);
end

dzdd = [];
correctionRHS = 0;
if opt.computeCorrection
    correctionRHS = 1;
dzdd = zeros(nx,1);
end

ax = []; 
if opt.computeCorrection
ax = cell(totalPredictionSteps,1);
end

if uSeedsProvided
    Ax = cell(totalPredictionSteps,1);
else 
Ax = cell(totalPredictionSteps,totalControlSteps);
end

av = [];
if opt.computeCorrection
av = cell(totalPredictionSteps,1);
end

if uSeedsProvided
    Av = cell(totalPredictionSteps,1);
else 
Av = cell(totalPredictionSteps,totalControlSteps);

end

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
        if ~uSeedsProvided && isempty(Ax{k-1,i})
            xRightSeeds = [{dzdd},Ax(k-1,1:i-1),{zeros(nx,nuSeed(i))}];
        elseif ~uSeedsProvided 
            xRightSeeds = [{dzdd},Ax(k-1,1:i)];
        elseif uSeedsProvided
            xRightSeeds = [{dzdd},Ax(k-1,1)];     
        end       
        
        seedSizes = cellfun(@(x)size(x,2),xRightSeeds);   
        xRightSeeds = cell2mat(xRightSeeds);
        
        if uSeedsProvided
            uRightSeeds = [zeros(uDims(i),correctionRHS),opt.uRightSeeds{i}];            
        else
            uRightSeeds = [zeros(uDims(i),correctionRHS+sum(nuSeed(1:i-1))),eye(uDims(i))];
        end

        xStart = x{k-1};
        
    elseif  k == 1
        
        xRightSeeds = zeros(nx,nuSeed(i));
        if uSeedsProvided
            uRightSeeds = opt.uRightSeeds{i};            
        else
            uRightSeeds = eye(uDims(1));

        end
        xStart = ss.state;
        
    else
        error('what?')
    end
    

    
    % Compute the Jacobian-vector products
    [xs{k},vs{k},Jac,convergence] = ss.step{k}(xStart,ui,'gradients',true,'xRightSeeds',xRightSeeds,'uRightSeeds',uRightSeeds,'simVars',simVars{k});
    
    converged(k) = convergence.converged;
    
    if ~givenRangeRHS
    xd{k} = (xs{k}-x{k});
    if withAlgs
        vd{k} = (vs{k}-v{k});
    end
    end
    
    % Extract the linearized model information
    if k >1
        [dzddAx] = mat2cell(Jac.xJ,nx,seedSizes);
        if uSeedsProvided
            [dzdd,Ax{k,1}] = deal(dzddAx{:});       
        else    
        [dzdd,Ax{k,1:i}] = deal(dzddAx{:});
        end
    else % k ==1
        [Ax{1,1}] = Jac.xJ;        
    end
    
    if opt.computeCorrection
    dzdd = dzdd -xd{k};
    
    ax{k}    = -dzdd;
    end
    
    if withAlgs
        if k >1
            nv = size(Jac.vJ,1);
            [dvddvsu] = mat2cell(Jac.vJ,nv,seedSizes);
            if uSeedsProvided
                [dvdd,Av{k,1}] = deal(dvddvsu{:});
            else
            [dvdd,Av{k,1:i}] = deal(dvddvsu{:});
            
            end
            if opt.computeCorrection
            av{k}    = vd{k}-dvdd;
            end
        else %k==1
            Av{1,1} = Jac.vJ;
            if opt.computeCorrection
            av{1}    = vd{k};
            end
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

