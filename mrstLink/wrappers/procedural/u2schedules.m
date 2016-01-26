function [ schedules ] = u2schedules( u,schedules,varargin)

opt = struct('uScale', [], 'fixedWells', []);
opt = merge_options(opt, varargin{:});

if isempty(opt.uScale)
    uScale = repmat({1},numel(u),1);
else
    uScale = opt.uScale;
end


schedulesSI = cellfun(@(ui,schedule,uScaleIt) controls2Schedule( ui,schedule,'uScale',uScaleIt, 'fixedWells', opt.fixedWells),...
               u,arrayfun(@(x){x}, schedules),uScale,...
               'UniformOutput',false);


schedulesSI = [schedulesSI{:}];

schedules = mergeSchedules(schedulesSI);





end

