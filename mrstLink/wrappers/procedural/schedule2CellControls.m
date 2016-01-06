function [ vals, type, control  ] = schedule2CellControls( schedule )
%
% extract schedule control information
%


useMrstSchedule = isfield(schedule.control(1), 'W');

if useMrstSchedule
    vals = arrayfun(@(c)zeros(numel(c.W),1),schedule.control,'UniformOutput',false);
    
    if nargout > 1
        control = arrayfun(@(c)cell(numel(c.W),1),schedule.control,'UniformOutput',false);
        type = arrayfun(@(c)cell(numel(c.W),1),schedule.control,'UniformOutput',false);
        
        for k = 1:numel(schedule.control)
            [ vals{k}, type{k}, control{k} ] = wells2Values(schedule.control(k).W);
        end
    else
        for k = 1:numel(schedule.control)
            [ vals{k}] = wells2Values(schedule.control(k).W);
        end
        
    end
else
    [ vals, type, control ] = getValuesFromSchedule(schedule) ;
end

%--------------------------------------------------------------------------
end




