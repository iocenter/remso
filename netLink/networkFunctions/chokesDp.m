function [obj] = chokesDp(forwardStates, schedule, p, netSol, nScale, pScale, varargin )
%CHOKESDP Calculates pressure drops of chokes in the network
    
    opt     = struct('ComputePartials',false, ...                                          
                     'leftSeed',[], ...
                     'activeComponents', struct('oil',1,'water',1,'gas',0,'polymer',0,'disgas',0,'vapoil',0,'T',0,'MI',0), ...
                     'fluid', [], ...
                    'dpFunction', @dpBeggsBrillJDJ, ...
                    'forwardGradient',true,...
                    'finiteDiff', true);
                 
                  
                     
    opt     = merge_options(opt, varargin{:});
    
    wellSols = cellfun(@(x)x.wellSol,forwardStates,'UniformOutput',false);    

    numSteps   = numel(forwardStates);
    
    obj = cell(1,numSteps);    
    lastStep = numSteps;
            
    wellSol = wellSols{lastStep};       
        
    netSol = runNetwork(netSol, wellSol, forwardStates{lastStep}, p, pScale, 'ComputePartials', opt.ComputePartials, 'dpFunction', opt.dpFunction, 'forwardGradient', opt.forwardGradient,'finiteDiff', opt.finiteDiff);   % running the network        
    
    
    obj{lastStep} = getChokesDp(netSol)./nScale; % returns pressure losses in chokes 

	if opt.ComputePartials && ~(size(opt.leftSeed,2)==0)
    	obj{lastStep}.jac = cellfun(@(x)opt.leftSeed*x,obj{lastStep}.jac,'UniformOutput',false);
	end
    
    if numSteps > 1
       obj0 =  obj{numSteps}*0;
      
        for step = 1:numSteps-1
            obj{step} = obj0;
        end
    end 
  
end

