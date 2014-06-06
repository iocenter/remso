function [ stateMrst ] = stateVector2stateMrst( stateVector,varargin)
%
%  write a state vector as a mrst state structure
%
%

opt = struct('xScale',[]);
opt = merge_options(opt, varargin{:});


% TODO: considering oil water
nx = numel(stateVector)/2;

if ~isempty(opt.xScale)
    stateVector = stateVector.*opt.xScale;
end


stateMrst.pressure = stateVector(1:nx);
stateMrst.s = [stateVector(nx+1:end),1-stateVector(nx+1:end)];


end