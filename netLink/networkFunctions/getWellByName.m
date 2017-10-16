function [ well ] = getWellByName(wellSol, wellName)
%getWellByName find a well by the name. Returns the wellSol mock object
%with the found well.
    i=1;
    while i<=numel(wellSol)
        if strcmp(wellSol(i).name, wellName)
            well = wellSol(i);
            return;
        end        
        i = i+1;
    end
    well = [];
end

