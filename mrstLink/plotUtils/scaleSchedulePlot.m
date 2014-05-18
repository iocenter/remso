function [v,schedulesSI] = scaleSchedulePlot(u,schedules,uScale,uScalePlot)

schedulesSI = cellfun(@(ui,schedule,uScaleIt) controls2Schedule( ui,schedule,'doScale',true,'uScale',uScaleIt),...
               u,arrayfun(@(x){x}, schedules),uScale,...
               'UniformOutput',false);
uSI = cellfun(@(x) schedule2Controls(x),schedulesSI,'UniformOutput',false);
schedulesPlot = cellfun(@(ui,schedulei,uScalei) controls2Schedule(ui,schedulei,'doScale',true,'uScale',uScalei),uSI,schedulesSI,uScalePlot,'UniformOutput',false);
[vals] =cellfun(@(sch)schedule2CellControls(sch),schedulesPlot,'UniformOutput',false);

% TODO: if more than cell in v?
v = cell2mat(cellfun(@(v)v{:},vals','UniformOutput',false));

schedulesSI = [schedulesSI{:}];





end