function [ netSol ] = createESPNetwork(ns, varargin)
% instantiate network inspired in an example used in the Master Thesis of
% Eirik Roysem (NTNU), but adapted to consider ESPs as shown in the paper
% 'Exploring the Potential of Model-Based Optimization in Oil Production 
% Gathering Networks with ESP-Produced, High Water Cut Wells'.
    opt = struct('withPumps', false);
    opt = merge_options(opt, varargin{:});

    nWells = length(ns.V);

    manifoldVert = newVertex(length(ns.V)+1, -1);
    ns = addVertex(ns, manifoldVert);    
    for i=1:nWells     
        isProducer = (ns.V(i).sign == -1);
        isInjector = (ns.V(i).sign == 1);        
        sign = isProducer*-1 + isInjector;
        
        if isProducer % production well infrastructure           
            inletBoosterVert = newVertex(length(ns.V)+1, sign);
            ns = addVertex(ns, inletBoosterVert);
            
            prodCasing =  newEdge(length(ns.E)+1,ns.V(i), inletBoosterVert, sign);
            prodCasing.pipeline = wellCasingSettings();
            prodCasing.stream = espStream();
            ns = addEdge(ns, prodCasing);
            
            if opt.withPumps
                outletBoosterVert = newVertex(length(ns.V)+1, sign);
                ns = addVertex(ns, outletBoosterVert);          
            
                booster = newEdge(length(ns.E)+1, inletBoosterVert, outletBoosterVert, sign);
                booster.stream = espStream();
                ns = addEdge(ns, booster, 'isEquipment', true);           
            end
            
            
            finalTubingVert = newVertex(length(ns.V)+1, sign);
            ns = addVertex(ns, finalTubingVert);
            
            if opt.withPumps
                prodTubing = newEdge(length(ns.E)+1,outletBoosterVert, finalTubingVert, sign);
            else
                prodTubing = newEdge(length(ns.E)+1,inletBoosterVert, finalTubingVert, sign);
            end
            prodTubing.units = 0; % METRIC=0, FIELD = 1,
            prodTubing.pipeline = wellTubingSettings();
            prodTubing.stream = espStream();
            ns = addEdge(ns, prodTubing); 
            
            if ~opt.withPumps %% if there is no pumps, include a choke in the well-head
                 outletChokeVert = newVertex(length(ns.V)+1, sign);
                 ns = addVertex(ns, outletChokeVert);                 
                 
                 choke = newEdge(length(ns.E)+1, finalTubingVert, outletChokeVert, sign);          
                 choke.stream = espStream();
                 ns = addEdge(ns, choke, 'isEquipment', true);
            end
                
            if opt.withPumps
                horizFlowline = newEdge(length(ns.E)+1, finalTubingVert, manifoldVert, sign);
            else
                horizFlowline = newEdge(length(ns.E)+1, outletChokeVert, manifoldVert, sign);
            end            
            horizFlowline.units = 0; % METRIC = 0, FIELD = 1           
            horizFlowline.pipeline = horizontalPipeSettings(ns.V(i).name);
            horizFlowline.stream = espStream();
            ns = addEdge(ns, horizFlowline);
            
        elseif isInjector % injection infrastruture
            
        end
    end
    
    inletSurfaceSepVert = newVertex(length(ns.V)+1, sign);    
    ns = addVertex(ns, inletSurfaceSepVert, 'isSink', true);
    
    flowlineRiser = newEdge(length(ns.E)+1, manifoldVert, inletSurfaceSepVert, 0);
    flowlineRiser.units = 0; % METRIC =0 , FIELD = 1
    flowlineRiser.pipeline = flowlineRiserSettings();
    flowlineRiser.stream = espStream();        
    ns = addEdge(ns, flowlineRiser);
    
    ns.boundaryCond = 5*barsa; % network boundary condition
                    
    netSol = ns;
end

function [str] = espStream() % default stream used in the example
    str = newStream('sg_gas', 0.65, ...  % air = 1
                    'oil_dens', 897,...  % kg/m^3
                    'water_dens', 1025.2, ... % kg/m^3  
                    'oil_visc', 0.00131, ... % Pa s
                    'water_visc', 0.00100); % Pa s
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
                  'len', 1000 , ... % in m
                  'ang', degtorad(0), ...% in rad
                  'temp',  convtemp(60,'C','K'));   % in C  
    elseif pipeOpt == 2        
        pipe = newPipeline('diam', 0.12, ... in %m
              'len', 1500 , ... % in m
              'ang', degtorad(5), ...% in rad
              'temp',  convtemp(60,'C','K'));   % in C 
    else        
        error('Standard pipeline should have been given !')
    end                  

end

function [pipe] = wellCasingSettings() %pipeW.dat
     pipe = newPipeline('diam', 152*milli*meter, ... in %m
                      'len', 213.3 , ... % in m
                      'ang', degtorad(90), ...  % in rad
                      'temp', convtemp(60,'C','K'));   % in K  

end

function [pipe] = wellTubingSettings() %pipeW.dat
     pipe = newPipeline('diam', 76*milli*meter, ... in %m
                      'len', 914.4 , ... % in m
                      'ang', degtorad(90), ...  % in rad
                      'temp', convtemp(60,'C','K'));   % in K  

end

function [pipe] = flowlineRiserSettings() %pipeR.dat
    pipe = newPipeline('diam', 0.24, ... in %m
                      'len', 2000 , ... % in m
                      'ang', degtorad(90), ...  % in rad
                      'temp',  convtemp(60,'C','K'));   % in K  
end
