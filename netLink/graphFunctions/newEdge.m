function edge = newEdge(nE, v1, v2, sign, stream, pipe)
   %% This functions creates a graph's edge representing 
   %% a flow line of the production network
   %% TODO: rewrite input list as optional parameters using opt and varargin
   
    streamGiven = (nargin == 5);
    if streamGiven
        str = stream;
    else
        str = defaultStream();
    end       
    
    pipeGiven = (nargin == 6);    
    if pipeGiven 
        pipeline = pipe;
    else
        pipeline = defaultPipeline();
    end
   
   
    %% edge of the production network which represents a flowline
    edge = struct(...
          'id',nE, ...
          'sign', sign, ...% -1 = 'UPHILL', 1 = 'DOWNHILL'         
          'name', strcat('e', num2str(nE)),... % name of the edge
          'vin', v1.id, ...   % start vertex
          'vout', v2.id, ...  % end vertex
          'units', 0, ...  % 0 = 'METRIC', 1 = 'FIELD'         
          'stream', str, ...
          'pipeline', pipeline, ...          
          'equipment', false, ...                              
          'integrationStep', 200*meter);  % integration step size            
end

function pipeline = defaultPipeline()
%% pipeline features
    pipeline = newPipeline('diam', 0.12, ... in %m
        'len', 1 , ... % in m
        'ang', deg2rad(90), ...  % in rad
        'temp', from_deg_C_to_deg_K(60)); % in K
        
end

function str = defaultStream()
    % black oil properties.
    str = struct('sg_gas', 0.83, ...      % gas specific gravity                                   
                 'oil_dens',0,...        % oil density
                 'water_dens', 1.02 ,...    % water density  
                 'oil_visc', 0, ...      % oil viscosity  
                 'water_visc', 0);    % water viscosity  
end
