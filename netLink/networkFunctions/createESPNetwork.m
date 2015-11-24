function [ netSol ] = createESPNetwork(ns)
% instantiate network inspired in an example used in the Master Thesis of
% Eirik Roysem (NTNU), but adapted to consider ESPs as shown in the paper
% 'Exploring the Potential of Model-Based Optimization in Oil Production 
% Gathering Networks with ESP-Produced, High Water Cut Wells'.
    nWells = length(ns.V);

    inletSubseaSepVert = newVertex(length(ns.V)+1, -1,-1);
    ns = addVertex(ns, inletSubseaSepVert);    
    for i=1:nWells     
        isProducer = (ns.V(i).sign == -1);
        isInjector = (ns.V(i).sign == 1);        
        sign = isProducer*-1 + isInjector;
        
        if isProducer % production well infrastructure      
            inletBoosterVert = newVertex(length(ns.V)+1, sign,sign);            
            ns = addVertex(ns, inletBoosterVert);
            
            prodTubing = newEdge(length(ns.E)+1, ns.V(i), inletBoosterVert, sign);
            prodTubing.units = 0; % METRIC=0, FIELD = 1,
            prodTubing.pipeline = wellRiserSettings();
            prodTubing.stream = eirikStream();
            ns = addEdge(ns, prodTubing);
            
            outletBoosterVert = newVertex(length(ns.V)+1, sign,sign);        
            ns = addVertex(ns, outletBoosterVert);
            
            booster = newEdge(length(ns.E)+1, inletBoosterVert, outletBoosterVert, sign);
            ns = addEdge(ns, booster, 'isPump', true);            
            
            horizFlowline = newEdge(length(ns.E)+1, outletBoosterVert, inletSubseaSepVert, sign);
            horizFlowline.units = 0; % METRIC = 0, FIELD = 1           
       

            horizFlowline.pipeline = horizontalPipeSettings(ns.V(i).name);
            horizFlowline.stream = eirikStream();
            ns = addEdge(ns, horizFlowline);
            
        elseif isInjector % injection infrastruture
            
        end
    end
    outletSubseaVert =  newVertex(length(ns.V)+1, sign, 0);
    ns = addVertex(ns, outletSubseaVert);
    
    subseaSeparator = newEdge(length(ns.E)+1, inletSubseaSepVert, outletSubseaVert, 0);
    ns = addEdge(ns, subseaSeparator, 'isSeparator', true);
    
    inletSurfaceSepVert = newVertex(length(ns.V)+1, sign, -1);
    inletSurfaceSepVert.pressure = 20; % in barsa 
    ns = addVertex(ns, inletSurfaceSepVert, 'isSink', true, 'isControllable', true);
    
    flowlineRiser = newEdge(length(ns.E)+1, outletSubseaVert, inletSurfaceSepVert, 0);
    flowlineRiser.units = 0; % METRIC =0 , FIELD = 1
    flowlineRiser.pipeline = flowlineRiserSettings();
    flowlineRiser.stream = eirikStream();        
    ns = addEdge(ns, flowlineRiser);
                    
    netSol = ns;
end

function [str] = eirikStream() % default stream used in the example
    str = newStream('sg_gas', 0.65, ...  % air = 1
                    'oil_dens', 700,...  % kg/m^3
                    'water_dens', 1000, ... % kg/m^3  
                    'oil_visc', 2, ... % cp
                    'water_visc', 0.9); % cp
end


function [pipe] = horizontalPipeSettings(wellName)
    if (strcmp(wellName, 'p1') || strcmp(wellName, 'p3') || ...
        strcmp(wellName, 'p4') || strcmp(wellName, 'p5'))    
        pipeOpt = 1;
    elseif strcmp(wellName, 'p2')
        pipeOpt = 2;
    else
        pipeOpt = -1;  
    end
    if pipeOpt == 1
        pipe = newPipeline('diam', 0.12, ... in %m
                  'len', 100 , ... % in m
                  'ang', degtorad(0), ...% in rad
                  'temp',  60);   % in C  
    elseif pipeOpt == 2        
        pipe = newPipeline('diam', 0.12, ... in %m
              'len', 150 , ... % in m
              'ang', degtorad(5), ...% in rad
              'temp',  60);   % in C 
    else        
        error('Standard pipeline should have been given !')
    end                  

end

function [pipe] = wellRiserSettings() %pipeW.dat
     pipe = newPipeline('diam', 0.12, ... in %m
                      'len', 1000 , ... % in m
                      'ang', degtorad(90), ...  % in rad
                      'temp',  60);   % in C   

end

function [pipe] = flowlineRiserSettings() %pipeR.dat
    pipe = newPipeline('diam', 0.24, ... in %m
                      'len', 2000 , ... % in m
                      'ang', degtorad(90), ...  % in rad
                      'temp',  60);   % in C   
end

