function [obj] = pumpsDp(forwardStates, schedule, p, netSol, nScale, numStages, pScale,  varargin )
%PumpsDp Calculates pressure drops of pumps in the network
    
    opt     = struct('ComputePartials',false, ...                                          
                     'leftSeed',[], ...                     
                     'dpFunction', @simpleDp);
                     
    opt     = merge_options(opt, varargin{:});
    
    wellSols = cellfun(@(x)x.wellSol,forwardStates,'UniformOutput',false);    

    numSteps   = numel(forwardStates);
    
    obj = cell(1,numSteps);
    
    step = numSteps;
        
    wellSol = wellSols{step};        
        
    netSol = runNetwork(netSol, wellSol, forwardStates{step}, p, pScale, 'ComputePartials', opt.ComputePartials, 'dpFunction', opt.dpFunction);   % running the network
   
    dpf = getChokesDp(netSol)./nScale; % dp in the pumps    
    
    obj{step} = dpf;

	if opt.ComputePartials && ~(size(opt.leftSeed,2)==0)
    	obj{step}.jac = cellfun(@(x)opt.leftSeed*x,obj{step}.jac,'UniformOutput',false);
	end
    
    if numSteps > 1
       obj0 =  obj{numSteps}*0;
      
        for step = 1:numSteps-1
            obj{step} = obj0;
        end
    end 
  
end

