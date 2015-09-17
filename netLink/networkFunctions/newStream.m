function str = newStream(varargin)
% this function creates a mock object with stream properties        
    if nargin ==0 
        str = defaultStream();
    else
        opt = struct('sg_gas', 0,'oil_dens',0, 'water_dens', 0 , 'oil_visc',0,'water_visc',0, 'gas_visc', 0);
        str = merge_options(opt, varargin{:});
    end
end

function str = defaultStream()
    % black oil properties.
    str = struct('sg_gas', 0.752, ...      % gas specific gravity                                   
                 'sg_oil', 0.7628, ...          % oil specific gravity
                 'gas_dens', 0.2, ...       % gas density
                 'oil_dens', 0.5,...        % oil density
                 'water_dens', 1.02 ,...    % water density  
                 'oil_visc', 1.3, ...      % oil viscosity  
                 'water_visc', 0.6, ...      % water viscosity  
                 'gas_visc', 0);            % gas viscosity
end

    