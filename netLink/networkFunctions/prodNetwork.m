function netSol = prodNetwork(wellSol, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%                                                                        %   
% Creates a simple production network given the installed production and %
% injection wells present in wellSol mock object                         %                                                                       %       
%                                                                        %   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    opt = struct('simpleNetwork',false, 'toyNetwork',false, 'eirikNetwork', false, 'espNetwork', false, 'satelliteWellsNetwork', false, 'withPumps', false);
    opt = merge_options(opt, varargin{:});

    netSol = initNetSolLocal(wellSol);       
    
    if opt.simpleNetwork
        netSol = createSimpleNetwork(netSol);  
    elseif opt.toyNetwork    
        netSol = createToyNetwork(netSol);
    elseif opt.eirikNetwork
        netSol = createEirikNetwork(netSol);
    elseif opt.espNetwork
        netSol = createESPNetwork(netSol, 'withPumps', opt.withPumps);
    elseif opt.satelliteWellsNetwork
        netSol = createSatelliteWellsNetwork(netSol);
    end
    
    % flow matrix
    netSol.M = createFlowMatrix(netSol);
    
    % init flow vectors
    netSol.qo = zeros(numel(netSol.E),1);
    netSol.qw = zeros(numel(netSol.E),1);
    netSol.qg = zeros(numel(netSol.E),1);        
    
    % init pressure vector
    netSol.pV = zeros(numel(netSol.V),1);
    netSol.pV(setdiff(netSol.Vsnk, netSol.VwInj)) = netSol.boundaryCond; % network boundary condition
end  

