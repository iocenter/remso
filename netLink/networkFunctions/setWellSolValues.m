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
    if ~opt.hasGas
        assert(norm(qGs)<100*eps)
        qGs = zeros(size(qGs));
    end
    
    pBHP = vertcat(wellSol.bhp);     
    
 
        
 
    if opt.ComputePartials  
        pressure  = forwardState.pressure;
        sW = forwardState.s(:,1);
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
        
        % intialize network with flows and pressures output from the reservoir simulation
        if opt.ComputePartials
            emptyADIp = repmat(pBHP(1)*0,numel(netSol.pV),1);
            emptyADIq = repmat(qOs(1)*0,numel(netSol.qo),1);
            
            netSol.pV = emptyADIp + netSol.pV;
            netSol.qo = emptyADIq + netSol.qo;
            netSol.qw = emptyADIq + netSol.qw;
            netSol.qg = emptyADIq + netSol.qg;
        end
    end
    netSol.pV(netSol.Vw) = pBHP;
    netSol.qo(netSol.VwProd) = qOs(netSol.VwProd);
    netSol.qw(netSol.VwProd) = qWs(netSol.VwProd);
    netSol.qg(netSol.VwProd) = qGs(netSol.VwProd);
    
        
end






