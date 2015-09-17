function edge = newEdge(nE, v1, v2, sign, stream, pipe)
   %% This functions creates a graph's edge representing 
   %% a flow line of the production network
   
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
          'vin', v1, ...   % start vertex
          'vout', v2, ...  % end vertex
          'units', 0, ...  % 0 = 'METRIC', 1 = 'FIELD'         
          'stream', str, ...
          'pipeline', pipeline, ...         
          'qoE', 0, ...
          'qwE', 0, ...
          'qgE', 0, ...
          'dpE', 0);
end

function pipeline = defaultPipeline()
    %% pipeline features
    pipeline =  struct('diameter', 0, ...      % inner pipe diameter in inches
                       'length', 0, ...        % total pipeline length (m to ft)
                       'angle', 0, ...         % inclination of pipe
                       'temperature', 0);      % average temperature in pipe
end



function str = defaultStream()
    % black oil properties.
    str = struct('sg_gas', 0.83, ...      % gas specific gravity                                   
                 'oil_dens',0,...        % oil density
                 'water_dens', 1.02 ,...    % water density  
                 'oil_visc', 0, ...      % oil viscosity  
                 'water_visc', 0);    % water viscosity  
end
