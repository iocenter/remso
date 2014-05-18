function [ newCons,nNew] = newConstraints(oldSet,newSet)
%
%  newCons are the constraints that are active in the new set, but are
%  not in the old set.  
%
%  nNew = number of new constraints
%

newCons = cellfun(@(q,w)~q&w,oldSet,newSet,'UniformOutput',false);

nNew = sum(cellfun(@sum,newCons));

end

