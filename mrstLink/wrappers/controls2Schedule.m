function [ schedule ] = controls2Schedule( u,schedule,varargin)
%
%  set the controls in the schedule.  Scale the variables
%
%

opt = struct('uScale',[]);
opt = merge_options(opt, varargin{:});


if ~isempty(opt.uScale)
    u = u.*opt.uScale;
end

vals = controls2CellControls(u,schedule);

schedule = cellControls2Schedule(vals,schedule);

