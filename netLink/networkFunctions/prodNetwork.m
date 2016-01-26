function netSol = prodNetwork(wellSol, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%                                                                        %   
% Creates a simple production network given the installed production and %
% injection wells present in wellSol mock object                         %                                                                       %       
%                                                                        %   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    opt = struct('simpleNetwork',false, 'toyNetwork',false, 'eirikNetwork', false, 'satelliteWellsNetwork', false);
    opt = merge_options(opt, varargin{:});

    netSol = initNetSolLocal(wellSol);       
    
    if opt.simpleNetwork
        netSol = createSimpleNetwork(netSol);  
    elseif opt.toyNetwork    
        netSol = createToyNetwork(netSol);
    elseif opt.eirikNetwork
        netSol = createEirikNetwork(netSol);
    elseif opt.satelliteWellsNetwork
        netSol = createSatelliteWellsNetwork(netSol);
    end
end  
