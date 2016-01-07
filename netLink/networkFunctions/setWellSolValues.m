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
    if opt.hasGas
        qGs = vertcat(wellSol.qGs);
    end
    pBHP = vertcat(wellSol.bhp);     
    
    pressure  = forwardState.pressure;
    sW = forwardState.s(:,1);
        
 
    if opt.ComputePartials
        if ~opt.hasGas
            [~, ~, qWs, qOs, pBHP, p] = initVariablesADI(pressure, sW, qWs, qOs, pBHP, p);
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
                p.jac{k} = pADI.jac{k};
                sW.jac{k} = sWADI.jac{k};
                rGH.jac{k} = rGHADI.jac{k};
            end
            
        end        
    end

    for i=1:length(wellSol)
        well =  getVertex(netSol, netSol.Vw(i));
        well.pressure =  pBHP(i);        
        well.qoV = qOs(i);
        
        if opt.hasGas        
                well.qgV = qGs(i);
        end
        well.qwV = qWs(i);    
        netSol = updateVertex(netSol, well);               
    end
    
    for j=1:length(netSol.Vc)
        vertControl = getVertex(netSol, netSol.Vc(j));
%         vertControl.pressure = p/barsa;
        vertControl.pressure = p*pScale;
        
        netSol = updateVertex(netSol, vertControl);
        
    end    
   
end

