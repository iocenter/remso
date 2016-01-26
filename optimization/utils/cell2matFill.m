function [ matrix ] = cell2matFill( cellMatrix,iDims,jDims)
% similar to cell2mat but fill the empty cells with zeros if the given size
%
%

if nargin < 2
    [di,dj] = size(cellMatrix);
    celldim = size(cellMatrix{1,1});
    iDims = celldim(1)*ones(di,1);
    jDims = celldim(2)*ones(1,dj);
end

if isrow(iDims)
    iDims = iDims';
end
cumI = cumsum([0;iDims(1:end-1)]);

if iscolumn(jDims)
    jDims = jDims';
end 
cumJ = cumsum([0,jDims(1:end-1)]);
   
cumI = repmat(num2cell(cumI),1,numel(jDims));
cumJ = repmat(num2cell(cumJ),numel(iDims),1);

[cellMatrixI,cellMatrixJ,cellMatrixV] = cellfun(@find,cellMatrix,'UniformOutput',false);

cellMatrixI = cellfun(@(ciM,cic)ciM+cic,cellMatrixI,cumI,'UniformOutput',false);
cellMatrixJ = cellfun(@(ciM,cic)ciM+cic,cellMatrixJ,cumJ,'UniformOutput',false);


i = vertcat(cellMatrixI{:});
j = vertcat(cellMatrixJ{:});
v = vertcat(cellMatrixV{:});

matrix = sparse(i,j,v,sum(iDims),sum(jDims));



end

