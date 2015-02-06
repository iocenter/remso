function [ target ] = concatenateMrstTargets(targets,sumLeftSeeds)

% creates a function target that is the concatenation of target1 and target2

%{

An mrst target takes an input the index j related to the step index, a
shootingSolN the solution at the end of the step, wellSol also at the end
of the step, and the schedule for the given step


%}

target = arroba(@vertcatFunctions,[1,2,3,4],{targets,sumLeftSeeds},true);


end

function [obj]  = vertcatFunctions(j,shootingSolN,wellSol,schedule,targets,sumLeftSeeds,varargin)

opt = struct('ComputePartials',false);
opt = merge_options(opt, varargin{:});



nt = numel(targets);

obj = cell(nt,1);

for k = 1:nt
    obj{k} = callArroba(targets(k),{j,shootingSolN,wellSol,schedule},...
        'ComputePartials',opt.ComputePartials);
end

if ~sumLeftSeeds
    obj = vertcat(obj{:});
else
    jac = obj{1}.jac;
    for k = 2:nt
        jac = cellfun(@(x,y)x+y,jac,obj{k}.jac,'UniformOutput',false);
    end
    obj = ADI(cell2mat(cellfun(@(x)double(x),obj,'UniformOutput',false)),jac);
end




end
