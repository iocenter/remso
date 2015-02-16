function [ uCells ] = schedules2CellControls(schedules,varargin)
%
%  extract the control information from the schedules
%
%
opt = struct('cellControlScales',[]);
opt = merge_options(opt, varargin{:});



n_sp = numel(schedules);
uCells = cell(n_sp,1);

for k = 1:n_sp
    
    if ~isempty(opt.cellControlScales)
        uCells{k} = schedule2Controls(schedules(k),'uScale',opt.cellControlScales{k});
    else
        uCells{k} = schedule2Controls(schedules(k));
    end
end




end

