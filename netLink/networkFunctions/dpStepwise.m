function [ dpTotal ] = dpStepwise(Eout,  qo, qw, qg, p, varargin)
%stepwiseDp calculates total pressure drop in a pipeline performing a 
% a forward or backward Euller integration

opt     = struct('dpFunction', @simpleDp, ...
                 'backward', false);

opt     = merge_options(opt, varargin{:});

dpTotal = 0*qw;
steps   = 0.*qw;

pipes = vertcat(Eout.pipeline);
pipeSizes = vertcat(pipes.len);
integrationSteps = vertcat(Eout.integrationStep);

condStop = steps >= pipeSizes;
while ~all(condStop)
    args = {Eout,  qo, qw, qg, p};
    dpF = arroba(opt.dpFunction, 1:numel(args) , [], true);
    dp = callArroba(dpF,args);    
    
    stepSize = integrationSteps;
    condStep = (steps + integrationSteps) > pipeSizes;
    if any(condStep)
        stepSize(condStep) = pipeSizes(condStep)-steps(condStep);
    end
    
    stepDp = dp.*stepSize;
    
    if opt.backward  % backward euller integration
        p = p + stepDp;
    else            % forward euller integration
        p = p - stepDp;
    end
    
    dpTotal = dpTotal + stepDp;
    steps = steps + stepSize;
    
    condStop = steps >= pipeSizes;    
end

end

