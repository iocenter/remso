function [obj] = chokesDp(forwardStates, schedule, netSol, nScale, p, varargin )
%CHOKESDP Calculates pressure drops of chokes in the network
    
    opt     = struct('ComputePartials',false, ...                                          
                     'leftSeed',[]);
                     
    opt     = merge_options(opt, varargin{:});
    
    wellSols = cellfun(@(x)x.wellSol,forwardStates,'UniformOutput',false);    

    numSteps   = numel(forwardStates);
    
    obj = cell(1,numSteps); %TODO: replace VwProd to Vw
    
    step = numSteps;
    
    
%     for step = 1:numSteps
        % pressure and saturaton vectors just used for place-holding
%         p  = forwardStates{step}.p;
%         sW = forwardStates{step}.s(:,1);
        
    wellSol = wellSols{step};
        
%         qWs  = vertcat(wellSol.qWs);
%         qOs  = vertcat(wellSol.qOs);        
%         pBHP = vertcat(wellSol.bhp);     
% 
%         if opt.ComputePartials
%             [~, ~, qWs, qOs, pBHP] = initVariablesADI(p, sW, qWs, qOs, pBHP);
%         end
        
    netSol = runNetwork(netSol, wellSol, forwardStates{step}, p, 'ComputePartials', opt.ComputePartials);   % running the network
    
    obj{step} = getChokesDp(netSol)./nScale; % returns pressure losses in chokes 

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

