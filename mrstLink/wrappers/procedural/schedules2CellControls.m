function [ uCells,Jacs ] = schedules2CellControls(schedules,varargin)
%
%  extract the control information from the schedules
%
%
opt = struct('cellControlScales',[]);
opt = merge_options(opt, varargin{:});


Jacs = [];
n_sp = numel(schedules);
if isempty(opt.cellControlScales)
    opt.cellControlScales = cell(n_sp,1);
end

if nargout > 1
    [uCells,Jacs] = arrayfun(@(s,i)schedule2Controls(s,'uScale',opt.cellControlScales{i}),schedules,(1:n_sp)','UniformOutput',false);
else
	[uCells] = arrayfun(@(s,i)schedule2Controls(s,'uScale',opt.cellControlScales{i}),schedules,(1:n_sp)','UniformOutput',false);
end



end

