function [ var ] = bringVariables( distVar,jobSchedule )
%
%  bring back variables to the client
%
%

assert(isa(distVar,'Composite'))


var = cell(size(jobSchedule.job2Work,1),1);


for w = 1:numel(distVar)
    distVarIt = distVar{w};
    var(jobSchedule.work2Job{w}) = distVarIt;
    
end


end
