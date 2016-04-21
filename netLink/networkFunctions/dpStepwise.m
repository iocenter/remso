function [ dpTotal ] = dpStepwise(Eout,  qo, qw, qg, p, varargin)
%stepwiseDp calculates total pressure drop in a pipeline performing a 
% a forward or backward Euller integration

opt     = struct('dpFunction', @simpleDp);

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
    p = p + stepDp;

    
    dpTotal = dpTotal + stepDp;
    steps = steps + stepSize;
    
    condStop = steps >= pipeSizes;    
end

end

