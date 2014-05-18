function [ vals, type, control  ] = schedule2CellControls( schedule )
%
% extract schedule control information
%


useMrstSchedule = isfield(schedule.control(1), 'W');

if useMrstSchedule
    vals = zeroVec(schedule);
    
    if nargout > 1
        control = zeroCell(schedule);
        type = zeroCell(schedule);
        
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



function zz = zeroVec(schedule)
zz = cell(1, numel(schedule.control));
for k = 1:numel(zz)
    zz{k} = zeros(numel(schedule.control(k).W), 1);
end

end

function zz = zeroCell(schedule)
zz = cell(1, numel(schedule.control));
for k = 1:numel(zz)
    zz{k} = cell( numel(schedule.control(k).W), 1);
end

end




