function [obj] = bhpDp(forwardStates, schedule, p, netSol, nScale, pScale, varargin )
%CHOKESDP Calculates pressure drops of chokes in the network
    
    opt     = struct('ComputePartials',false, ...                                          
                     'leftSeed',[]);
                     
    opt     = merge_options(opt, varargin{:});
    
    wellSols = cellfun(@(x)x.wellSol,forwardStates,'UniformOutput',false);    

    numSteps   = numel(forwardStates);
    
    obj = cell(1,numSteps);
    
    step = numSteps;       
    
    wellSol = wellSols{step};
        
    netSol = runNetwork(netSol, wellSol, forwardStates{step}, p, pScale, 'ComputePartials', opt.ComputePartials);   % running the network
    
    obj{step} = getBhpResidual(wellSol, netSol)./nScale; % returns bhp pressure residual

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

