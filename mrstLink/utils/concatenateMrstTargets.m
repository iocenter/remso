function [ target ] = concatenateMrstTargets(targets,sumLeftSeeds,outDimens)

% creates a function target that is the concatenation of target1 and target2

%{

An mrst target takes as input the pair (forwardStates,schedule) obtained from
simulation

sumLeftSeeds must be set to true only if 'targets' include their
corresponding leftSeeds at construction.

%}

if nargin <3
    outDimens = 0;
end

target = arroba(@vertcatFunctions,[1,2,3],{targets,sumLeftSeeds,outDimens},true);


end

function [obj]  = vertcatFunctions(forwardStates,schedule,p,targets,sumLeftSeeds,outDimens,varargin)

opt = struct('ComputePartials',false,'leftSeed',[]);
opt = merge_options(opt, varargin{:});

if size(opt.leftSeed,2) > 0
    assert(~sumLeftSeeds, 'The left seeds were given before!');
    sumLeftSeeds = true;
    leftSeed = mat2cell(opt.leftSeed,size(opt.leftSeed,1),outDimens);
    targets = num2cell(targets);
    obj = cellfun(@(ti,lSi)callArroba(ti,{forwardStates,schedule,p},'ComputePartials',opt.ComputePartials,'leftSeed',lSi),targets,leftSeed,'UniformOutput',false);
else
    obj = arrayfun(@(ti)callArroba(ti,{forwardStates,schedule,p},'ComputePartials',opt.ComputePartials),targets,'UniformOutput',false);
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


