function [ stateVector,Jac ] = stateMrst2stateVector( stateMrst,varargin )
%
%
%  flatten a mrst state as a state vector
%  TODO: reimplement using finalStepVars!
%
%
opt = struct('xScale',[],...
    'activeComponents',struct('oil',1,'water',1,'gas',0,'polymer',0,'disgas',0,'vapoil',0,'T',0,'MI',0),...% default OW
    'fluid',[],...
    'system',[],...
    'partials',false);

opt = merge_options(opt, varargin{:});

comp = opt.activeComponents;


Jac = [];
if ~comp.gas && ~comp.polymer && ~(comp.T || comp.MI)
    
    stateVector = [stateMrst.pressure;stateMrst.s(:,1)];
    
    
    if opt.partials
        Jac = speye(numel(stateVector));
    end
    
elseif comp.gas && comp.oil && comp.water
    
    % The transformation function may be given as an input and
    % generalized    
    [ p,sW,rGH ] = stateMrst2statePsWrGH(stateMrst,opt.fluid,opt.system,'partials',opt.partials);
    
    stateVector = [double(p);
                   double(sW);
                   double(rGH)];
               
    if opt.partials
        obj = [p;sW;rGH];
        obj = cat(obj);
        Jac = obj.jac{1};
    end

else
    error('Not implemented for current activeComponents');
end



if ~isempty(opt.xScale)
    stateVector = stateVector./opt.xScale;
    if opt.partials
        Jac = bsxfun(@ldivide,opt.xScale,Jac);
    end
end


end