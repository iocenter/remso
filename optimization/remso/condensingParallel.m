function varargout = condensingParallel(x,u,v,ss,jobSchedule,simVars)
% Apply the Lift-Opt trick and get the correction and predictor matrices.
%
% SYNOPSIS:
%  [xd,vd,ax,Ax,av,Av] = condensing(x,u,v,ss)
%  [xd,vd,ax,Ax,av,Av] = condensing(x,u,v,ss, 'pn', pv, ...)
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


% Pepare inputs

label = 'condensing';
fprintf(label);

totalPredictionSteps = getTotalPredictionSteps(ss);
totalControlSteps = numel(u);

nx = numel(ss.state);
nv = ss.nv;
nu = numel(u{1});

withAlgs = nv>0;

xd = cell(totalPredictionSteps,1);
vd = cell(totalPredictionSteps,1);

ax = cell(totalPredictionSteps,1);
Ax = cell(totalPredictionSteps,totalControlSteps);
av = cell(totalPredictionSteps,1);
Av = cell(totalPredictionSteps,totalControlSteps);


if isempty(simVars)
    error('run simulateSystem first!')
elseif iscell(simVars)
    if isempty(simVars{1})
        error('run simulateSystem first!')
    end
else
    simVars = bringVariables(simVars,ss.jobSchedule);
end
workerCondensingSchedule = jobSchedule.workerCondensingSchedule;
uStart = jobSchedule.uStart;
step = ss.step;
ci = ss.ci;
state0 = ss.state;

spmd
    
    computeCorrection = min(workerCondensingSchedule) == 0;
    computeControl =  false;
    
    if computeCorrection
        dzdd = zeros(nx,1);
    else
        dzdd = [];
    end
    
    xdW = cell(totalPredictionSteps,1);
    vdW = cell(totalPredictionSteps,1);
    
    axW = cell(totalPredictionSteps,1);
    AxW = cell(totalPredictionSteps,totalControlSteps);
    avW = cell(totalPredictionSteps,1);
    AvW = cell(totalPredictionSteps,totalControlSteps);
    
    if computeCorrection
        kStart  = 1;
    else
        kStart = min(uStart(workerCondensingSchedule));
    end
    
    for k = kStart:totalPredictionSteps
        
        
        i = callArroba(ci,{k});
        ui = u{i};
        
        % Prepare the RHS for the condensing technique
        gradients = true;
        if k >1
            
            if isempty(AxW{k-1,i}) && any(ismember(workerCondensingSchedule,i))
                xRightSeeds = [{dzdd},AxW(k-1,1:i-1),{zeros(nx,nu)}];
                computeControl = true;
            else
                xRightSeeds = [{dzdd},AxW(k-1,1:i)];
            end
            
            seedSizes = fSeedSizes(xRightSeeds);
            
            uRightSeeds = fuRightSeedsAlloc(seedSizes,nu);
            
            
            if any(ismember(workerCondensingSchedule,i))
                uRightSeeds{i+1} = eye(nu);
            end
            
            xRightSeeds = cell2mat(xRightSeeds);
            uRightSeeds = cell2mat(uRightSeeds);
            
            xStart = x{k-1};
            
        elseif  k == 1
            computeControl = any(ismember(workerCondensingSchedule,i));
            
            if computeControl
                xRightSeeds = zeros(nx,nu);
                uRightSeeds = eye(nu);
            else % then computeCorrection
                xRightSeeds = zeros(nx,0);
                uRightSeeds = zeros(nu,0);
                gradients = false;
            end
            
            xStart = state0;
            
        else
            error('what?')
        end
        
        
        
        % Compute the Jacobian-vector products
        [xs,vs,Jac] = callArroba(step{k},{xStart,ui},...
            'gradients',gradients,...
            'xRightSeeds',xRightSeeds,...
            'uRightSeeds',uRightSeeds,...
            'simVars',simVars{k});
        
        
        if computeCorrection
            xdW{k} = (xs-x{k});
            if withAlgs
                vdW{k} = (vs-v{k});
            end
        end
        % Extract the linearized model information
        if k >1
            [dzddAx] = mat2cell(Jac.xJ,nx,seedSizes);
            [dzdd,AxW{k,1:i}] = deal(dzddAx{:});
        else % k ==1
            if computeControl
                [AxW{1,1}] = Jac.xJ;
            end
            if computeCorrection
                dzdd = zeros(nx,1);
            end
        end
        
        if computeCorrection
            dzdd = dzdd -xdW{k};
            axW{k}    = -dzdd;
        end
        
        if withAlgs
            if k >1
                [dvddvsu] = mat2cell(Jac.vJ,nv,seedSizes);
                [dvdd,AvW{k,1:i}] = deal(dvddvsu{:});
                
                if computeCorrection
                    avW{k}    = vdW{k}-dvdd;
                end
            else %k==1
                if computeControl
                    AvW{1,1} = Jac.vJ;
                end
                if computeCorrection
                    avW{1}    = vdW{k};
                end
            end
            
        end
        
        
    end
    
end
cs = jobSchedule.clientCondensingSchedule;

 varargout{1} = xdW{cs.correction};
 varargout{3} = axW{cs.correction};
 
 if withAlgs
    varargout{2} = vdW{cs.correction};
    varargout{5} = avW{cs.correction};
 end
 
 for k = 1:totalControlSteps
    Axk = AxW{cs.control(k)};
    Ax(:,k) = Axk(:,k);
   
    if withAlgs
        Avk = AvW{cs.control(k)};
        Av(:,k) = Avk(:,k);
    end
 end
 
  varargout{4} = Ax;
  varargout{6} = Av;
  
  fprintf(repmat('\b', 1, numel(label)));

end

function seedSizes = fSeedSizes(xRightSeeds)
    seedSizes = cellfun(@(x)size(x,2),xRightSeeds);
end
function uRightSeeds = fuRightSeedsAlloc(seedSizes,nu)
    uRightSeeds = arrayfun(@(z)zeros(nu,z),seedSizes,'UniformOutput',false);
end
            
            
            
