function [ set ] = unionActiveSets(set1,set2)
%  set = set1 or set2 
%

set = cellfun(@(q,w)or(w,q),set1,set2,'UniformOutput',false);


end

