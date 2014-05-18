function [ u ] = schedule2Controls(schedule,varargin)
%
% extract the schedule control information and scale it
%
%

opt = struct('doScale',false,'uScale',[]);
opt = merge_options(opt, varargin{:});


[ vals ] = schedule2CellControls(schedule);
u = cellControls2Controls(vals);

if opt.doScale
    if isempty(opt.uScale)
        [ scheduleScale ] = schedulesScaling( schedule );
        [ uScale ] = schedule2Controls(scheduleScale);
        u = u./uScale;
    else
        u = u./opt.uScale;
    end
end


end