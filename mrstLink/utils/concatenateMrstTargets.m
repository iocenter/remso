function [ target ] = concatenateMrstTargets(targets,sumLeftSeeds)

% creates a function target that is the concatenation of target1 and target2

%{

An mrst target takes an input the index j related to the step index, a
shootingSolN the solution at the end of the step, wellSol also at the end
of the step, and the schedule for the given step


%}

target = arroba(@vertcatFunctions,[1,2],{targets,sumLeftSeeds},true);


end

function [obj]  = vertcatFunctions(forwardStates,schedule,targets,sumLeftSeeds,varargin)

opt = struct('ComputePartials',false);
opt = merge_options(opt, varargin{:});



nt = numel(targets);

obj = cell(nt,1);

for k = 1:nt
    obj{k} = callArroba(targets(k),{forwardStates,schedule},...
        'ComputePartials',opt.ComputePartials);
end

if ~sumLeftSeeds
    obj = cellfun(@(varargin)vertcat(varargin{:}),obj{:},'UniformOutput',false);
else
    
    val = cellfun(@(t) extractField('val',t),obj,'UniformOutput',false);
    val = cellfun(@(varargin)vertcat(varargin{:}),val{:},'UniformOutput',false);
    
    jac = cellfun(@(t) extractField('jac',t),obj,'UniformOutput',false);
	jac = cellfun(@(varargin)sumCells(varargin{:}),jac{:},'UniformOutput',false);

   
    obj = cellfun(@(v,j)ADI(v,j),val,jac,'UniformOutput',false);
end




end

function out = extractField(field,t)
    out = cellfun(@(ti)ti.(field),t,'UniformOutput',false);
end

function out = sumCells(varargin)
    args = varargin;
    out = cellfun(@(varargin)sum2(varargin{:}),args{:},'UniformOutput',false);
end


function out = sum2(varargin)
%    out = sum(cat(3,varargin{:}),3);  this doesn't work for sparse
%    matrices

    matrixDim = size(varargin{1});
    nM = numel(varargin);
    
    S = repmat(speye(matrixDim(1)),1,nM);
    M = vertcat(varargin{:});

    out = S*M;

end


