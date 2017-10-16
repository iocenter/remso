function str = newStream(varargin)
% this function creates a mock object with stream properties        
    if nargin ==0 
        str = defaultStream();
    else
        opt = defaultStream();
        str = merge_options(opt, varargin{:});
    end
end

function str = defaultStream()
    % black oil properties.
    str = struct('sg_gas', 0.752, ...      % gas specific gravity                                   
                 'sg_oil', 0.7628, ...          % oil specific gravity
                 'gas_dens', 0.8, ...       % gas density
                 'oil_dens', 800,...        % oil density
                 'water_dens', 1000 ,...    % water density  
                 'oil_visc', 1.3, ...      % oil viscosity  
                 'water_visc', 0.6, ...      % water viscosity  
                 'gas_visc', 0);            % gas viscosity
end

    