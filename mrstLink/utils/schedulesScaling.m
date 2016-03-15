function [ schedules ] = schedulesScaling( schedules,varargin )
%
% fill the schedules mock object with the scaling information provided in
% the options
%
%

opt = struct('RATE',meter^3/day,'GRAT',meter^3/day,'ORAT',meter^3/day,'WRAT',meter^3/day,'LRAT',meter^3/day,'RESV',0,'BHP',barsa, 'FREQ', 1);
opt     = merge_options(opt, varargin{:});


n_sp = numel(schedules);


for k = 1:n_sp

    [ vals,type,control] = schedule2CellControls(schedules(k));
   
    for j = 1:numel(vals)  % number of set of well controls
        for i = 1:numel(vals{j}) % number of wells
            switch type{j}{i}
                case {'inj',1}
                    vals{j}(i) = opt.(control{j}{i});
                case {'prod',-1}
                    vals{j}(i) = opt.(control{j}{i});
                otherwise
                    error(['Cannot handle well type: ', type{j}]);
            end
        end
    end
    schedules(k) = cellControls2Schedule(vals,schedules(k));
    
end



end

