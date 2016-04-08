function [netSol] = setWellSolValues(netSol, wellSol, forwardState, p, pScale, varargin)
%SETWELLSOLVALUES set wellSol values in the network
%TODO: handle gas phase flow.

    opt     = struct('ComputePartials',false,...
                     'activeComponents',struct('oil',1,'water',1,'gas',0,'polymer',0,'disgas',0,'vapoil',0,'T',0,'MI',0), ...
                     'hasGas', false, ...
                     'fluid',[]);                     
    opt     = merge_options(opt, varargin{:});
    
    qWs  = vertcat(wellSol.qWs);
    qOs  = vertcat(wellSol.qOs);        
    qGs = vertcat(wellSol.qGs);
    
    pBHP = vertcat(wellSol.bhp);     
    
    pressure  = forwardState.pressure;
    sW = forwardState.s(:,1);
        
 
    if opt.ComputePartials        
        if ~opt.hasGas
            [~, ~, qWs, qOs, pBHP, p] = initVariablesADI(pressure, sW, qWs, qOs, pBHP, p);                       
            % initializing empty ADI for the gas phase
            qGs = 0*qWs;
        elseif opt.hasGas
            
            % The transformation function may be given as an input and
            % generalized
            disgas = opt.activeComponents.disgas;
            vapoil = opt.activeComponents.vapoil;
            fluid  =  opt.fluid;
            finalState = forwardState;
            
            [ press,sW,rGH ] = stateMrst2statePsWrGH(finalState,fluid,disgas,vapoil,'partials',opt.ComputePartials);
                     
            
            % instantiating jacobians with right values and dimensions.
            % Observe that the independet variables are p,sw,x,qw,qo,qg,bhp
            % and not what the line below suggest!
            [pADI,sWADI,rGHADI,qWs,qOs,qGs,pBHP, ~] = initVariablesADI(double(press),double(sW),double(rGH),qWs,qOs,qGs,pBHP, p);
            
            % revert jacobian given by stateMrst2statePsWrGH
            for k = 4:numel(pADI.jac)
                press.jac{k} = pADI.jac{k};
                sW.jac{k} = sWADI.jac{k};
                rGH.jac{k} = rGHADI.jac{k};
            end
            
        end
        %% initializing constants of network boundary conditions with empty jacobian
        surfaceSinks = setdiff(netSol.Vsnk,netSol.VwInj);
        for k=1:numel(surfaceSinks)
            vertSink = getVertex(netSol, surfaceSinks(k));
            pressureVal = vertSink.pressure;
            vertSink.pressure = qOs(1)*0 + pressureVal;
            netSol = updateVertex(netSol, vertSink);
        end
    end

    % update well vertices with corresponding flows and pressures
    for i=1:length(wellSol)
        well =  getVertex(netSol, netSol.Vw(i));
        well.pressure =  pBHP(i);        
        well.qoV = qOs(i);        
        well.qgV = qGs(i);        
        well.qwV = qWs(i);    
        netSol = updateVertex(netSol, well);               
    end
    
    for j=1:numel(netSol.Vc) % controllable vertices
        vertControl = getVertex(netSol, netSol.Vc(j));
        vertControl.pressure = p*pScale;   %% TODO: generalize this using the field 'control' the vertex mock object        
        
        netSol = updateVertex(netSol, vertControl);
        
    end 
    
    %%TODO: update set of controllable edges in createESPNetwork
%     for k=1:numel(netSol.Ec) % controllable edges        
%        edgeControl = getEdge(netSol, netSol.Ec(k));
%        edgeControl.control = p(k);  
%        
%        netSol = updateEdge(netSol, edgeControl);
%     end

   
   
end

