function [ stateMrst ] = stateVector2stateMrst( stateVector,varargin)
%
%  write a state vector as a mrst state structure
%
%

opt = struct('doScale',false,'xScale',[]);
opt = merge_options(opt, varargin{:});


% TODO: considering oil water
nx = numel(stateVector)/2;


if opt.doScale
    if isempty(opt.xScale)
        [stateScale] = stateScaling(stateVector);
        [ xScale ] = stateVector2stateMrst( stateScale );
        stateVector = stateVector.*xScale;
    else
        stateVector = stateVector.*opt.xScale;
    end
end



stateMrst.pressure = stateVector(1:nx);
stateMrst.s = [stateVector(nx+1:end),1-stateVector(nx+1:end)];


end