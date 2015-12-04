function [obj] = pumpsDp(forwardStates, schedule, p, netSol, nScale, numStages,  varargin )
%CHOKESDP Calculates pressure drops of pumps in the network
    
    opt     = struct('ComputePartials',false, ...                                          
                     'leftSeed',[]);
                     
    opt     = merge_options(opt, varargin{:});
    
    wellSols = cellfun(@(x)x.wellSol,forwardStates,'UniformOutput',false);    

    numSteps   = numel(forwardStates);
    
    obj = cell(1,numSteps);
    
    step = numSteps;
        
    wellSol = wellSols{step};        
        
    netSol = runNetwork(netSol, wellSol, forwardStates{step}, p, [], 'ComputePartials', opt.ComputePartials);   % running the network
    
    vw = getVertex(netSol, netSol.VwProd);
    ew = getEdge(netSol, vertcat(vw.Eout));
    
    qf = vertcat(ew.qoE) + vertcat(ew.qwE);  % flows in the pumps    
   
    dpf = getChokesDp(netSol); % dp in the pumps    
   
    inletStr = vertcat(ew.stream);
    mixtureDen = (vertcat(inletStr.oil_dens) + vertcat(inletStr.water_dens))./2;  % density of the mixture
    
    dhf= pump_dh(dpf, mixtureDen); % dh in the pumps
    
    [freq] = pump_eq_system_explicit(qf, dhf, 60, numStages);  % solves a system of equations to obtain frequency, flow and dh at 60Hz
    
    obj{step} = freq;
    
%      obj{step} = getChokesDp(netSol)./nScale; % returns pressure losses in chokes 

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

