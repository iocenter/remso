function [ stateVector ] = stateMrst2stateVector( stateMrst,varargin )
%
%
%  flatten a mrst state as a state vector
%
%
opt = struct('doScale',false,'xScale',[]);
opt = merge_options(opt, varargin{:});


stateVector = [stateMrst.pressure;stateMrst.s(:,1)]; % TODO: only for oil-Water systems


if opt.doScale
    if isempty(opt.xScale)
        [stateScale] = stateScaling(stateMrst);
        [ xScale ] = stateVector2stateMrst( stateScale );
        stateVector = stateVector./xScale;
    else
        stateVector = stateVector./opt.xScale;
    end
end


end