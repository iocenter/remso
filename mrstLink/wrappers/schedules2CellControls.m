function [ uCells ] = schedules2CellControls(schedules,varargin)
%
%  extract the control information from the schedules
%
%
opt = struct('doScale',false,'cellControlScales',[]);
opt = merge_options(opt, varargin{:});



n_sp = numel(schedules);
uCells = cell(n_sp,1);

for k = 1:n_sp
    
    if opt.doScale
        uCells{k} = schedule2Controls(schedules(k),'doScale',true,'uScale',opt.cellControlScales{k});
    else
        uCells{k} = schedule2Controls(schedules(k),'doScale',false);
    end
end




end

