function [ schedule ] = wellSol2Schedule( wellSol,schedule )

% poulate the well control values with the values given in the wellSol

[ nW ] = getnW( schedule );

if isfield(schedule.control,'W')

    for w = 1:nW
        schedule.control.W(w).val = wellSol(w).val;
    end
else
    error('not implemented');
end

end

