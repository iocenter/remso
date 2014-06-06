function [ u ] = schedule2Controls(schedule,varargin)
%
% extract the schedule control information and scale it
%
%

opt = struct('uScale',[]);
opt = merge_options(opt, varargin{:});


[ vals ] = schedule2CellControls(schedule);
u = cellControls2Controls(vals);

if ~isempty(opt.uScale)
    u = u./opt.uScale;
end



end