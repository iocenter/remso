function [ stateVector ] = stateMrst2stateVector( stateMrst,varargin )
%
%
%  flatten a mrst state as a state vector
%
%
opt = struct('xScale',[]);
opt = merge_options(opt, varargin{:});


stateVector = [stateMrst.pressure;stateMrst.s(:,1)]; % TODO: only for oil-Water systems

if ~isempty(opt.xScale)
    stateVector = stateVector./opt.xScale;
end


end