function [ stateMrst,Jac ] = stateVector2stateMrst( stateVector,varargin)
%
%  write a state vector as a mrst state structure
%
%

opt = struct('xScale',[],...
    'activeComponents',struct('oil',1,'water',1,'gas',0,'polymer',0,'disgas',0,'vapoil',0,'T',0,'MI',0),...% default OW
    'fluid',[],...
    'partials',false);

opt = merge_options(opt, varargin{:});

comp = opt.activeComponents;

if ~isempty(opt.xScale)
    stateVector = stateVector.*opt.xScale;
end

if ~comp.gas && ~comp.polymer && ~(comp.T || comp.MI)
    
    nx = numel(stateVector)/2;
    
    stateMrst.pressure = stateVector(1:nx);
    stateMrst.s = [stateVector(nx+1:end),1-stateVector(nx+1:end)];
    
    if opt.partials
        if ~isempty(opt.xScale)
            Jac = bsxfun(@times,speye(2*nx),opt.xScale');
        else
            Jac = speye(2*nx);
        end
    else
        Jac = [];
    end
    
    
    
elseif comp.gas && comp.oil && comp.water
    
    nx = numel(stateVector)/3;
    
    
    p = stateVector(1:nx);
    sW = stateVector(nx+1:2*nx);
    rGH = stateVector(2*nx+1:end);
    
    % The transformation function may be given as an input and
    % generalized
    disgas = opt.activeComponents.disgas;
    vapoil = opt.activeComponents.vapoil;
    [ stateMrst,Jac ] = statePsWrGH2stateMRST( p,sW,rGH,opt.fluid,disgas,vapoil,'partials',opt.partials);
    
    if opt.partials
        Jac = cat(2,Jac{:});
        if ~isempty(opt.xScale)
            Jac = bsxfun(@times,Jac,opt.xScale');
        end
    end
    
else
    error('Not implemented for current activeComponents');
end


end