function [x] = setStateValues(stateValue,varargin)
%
%  Given a structure with physical state values, i.e., (pressure,sW), provide a remso
%  state with these values.  If a optional 'x' value is given, the state
%  will be modified according to the optional parameters


opt = struct('x',[],'cells',[],'nCells',[],'xScale',[]);
opt = merge_options(opt, varargin{:});


if isempty(opt.x)  %%
    if numel(stateValue.pressure) == 1
        
        x = [ones(opt.nCells,1)*stateValue.pressure;
            ones(opt.nCells,1)*stateValue.sW];
        
    else %% assuming that correct dimensions are given!
        
        x = [stateValue.pressure;
            stateValue.sW];  
    end
else
    
    if ~isempty(opt.xScale)
        x = opt.x.*opt.xScale;
    else
        x = opt.x;
    end
    
    if  isfield(stateValue,'pressure')
        x(opt.cells) = stateValue.pressure;
    end
    if  isfield(stateValue,'sW')
        x(opt.cells+opt.nCells) = stateValue.sW;
    end
    
end


if ~isempty(opt.xScale)
    x = x./opt.xScale;
end



end
