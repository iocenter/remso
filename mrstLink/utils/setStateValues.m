function [x] = setStateValues(stateValue,varargin)
%
%  Given a structure with physical state values, i.e., (pressure,sW) or
%  (pressure,sW,rGH) for 2-phase or 3-phase respectively, provide a remso
%  state with these values.  If a optional 'x' value is given, the state
%  will be modified according to the optional parameters

%% TODO: This function is handy must should be deprecated


opt = struct('x',[],'cells',[],'nCells',[],'xScale',[],'scaling',[]);
opt = merge_options(opt, varargin{:});

if ~isempty(opt.xScale) && ~isempty(opt.scaling)
    error('use only one scaling method');
end
if isempty(opt.xScale)
    opt.xScale = 1;
end
if isempty(opt.scaling)
    opt.scaling = struct ('p',1,'s',1,'rGH',1);
end


if isempty(opt.x)  %%
    if numel(stateValue.pressure) == 1
        if isfield(stateValue,'rGH')
            x = [ones(opt.nCells,1)*stateValue.pressure;
                ones(opt.nCells,1)*stateValue.sW;
                ones(opt.nCells,1)*stateValue.rGH];
        else
            x = [ones(opt.nCells,1)*stateValue.pressure;
                ones(opt.nCells,1)*stateValue.sW];
        end
        
    else %% assuming that correct dimensions are given!
        if isfield(stateValue,'rGH')
            x = [stateValue.pressure;
                stateValue.sW;
                stateValue.rGH];
        else
            x = [stateValue.pressure;
                stateValue.sW];
        end
        
    end
else
    
    x = opt.x.*opt.xScale;
    if isfield(stateValue,'rGH')
        opt.nCells = numel(x)/3;
        x = x.*[ones(opt.nCells,1)*opt.scaling.p;
            ones(opt.nCells,1)*opt.scaling.s;
            ones(opt.nCells,1)*opt.scaling.rGH];
    else
        opt.nCells = numel(x)/2;
        x = x.*[ones(opt.nCells,1)*opt.scaling.p;
            ones(opt.nCells,1)*opt.scaling.s];
    end
    
    if  isfield(stateValue,'pressure')
        x(opt.cells) = stateValue.pressure;
    end
    if  isfield(stateValue,'sW')
        x(opt.cells+opt.nCells) = stateValue.sW;
    end
    if isfield(stateValue,'rGH')
        x(opt.cells+2*opt.nCells) = stateValue.rGH;
    end
    
end


x = x./opt.xScale;
if isfield(stateValue,'rGH')
    x = x./[ones(opt.nCells,1)*opt.scaling.p;
        ones(opt.nCells,1)*opt.scaling.s;
        ones(opt.nCells,1)*opt.scaling.rGH];
else
    x = x./[ones(opt.nCells,1)*opt.scaling.p;
        ones(opt.nCells,1)*opt.scaling.s];
end



end
