function [ violation,lowActive,upActive ] = checkConstraintViolation(x,lb,lowActive,ub,upActive,varargin)
%  violation = sum of all the constraints violation (norm1 of the violation)

opt = struct('first',0);
opt = merge_options(opt, varargin{:});



if opt.first == 0
    geqZeroNorm = @(x)max([x;0]);
    
    % if you want infinity norm
    violation = max([cellfun(@(var,upBound,act) geqZeroNorm(var(act) -upBound(act)),x,ub,upActive);...
        cellfun(@(var,lowBound,act)geqZeroNorm(lowBound(act)-var(act)),x,lb,lowActive)]);
    %if you want norm_1
    %violation = sum([cellfun(@(var,upBound,act) geqZeroNorm(var(act) -upBound(act)),x,ub,upActive);...
    %                 cellfun(@(var,lowBound,act)geqZeroNorm(lowBound(act)-var(act)),x,lb,lowActive)]);
    
else
    
    geqZero= @(x)max(x,0);
    
    violationUp  = cellfun(@(var,upBound,act) geqZero(var(act) -upBound(act)),x,ub,upActive,'UniformOutput',false);
    vuN = numel(violationUp);
    
    violationLow = cellfun(@(var,lowBound,act)geqZero(lowBound(act)-var(act)),x,lb,lowActive,'UniformOutput',false);
    
    violation = [violationUp;violationLow];
    vioDim = cellfun(@numel,violation);
    vioSum = sum(vioDim);
    vioAct = min([vioSum,opt.first]);
    
    violation = cell2mat(violation);
    
    [val,ord] = sort(violation,1,'descend');
    
    
    active = false(vioSum,1);
    active(ord(1:vioAct)) = true;
    
    
    active = mat2cell(active,vioDim,1);
    
    activeU = active(1:vuN);
    activeL = active(vuN+1:end);
    
    
    for k=1:numel(upActive)
        upActive{k}(upActive{k}) = activeU{k};
    end
    for k=1:numel(lowActive)
        lowActive{k}(lowActive{k}) = activeL{k};
    end
    if ~isempty(violation)
        violation = violation(1);
    else
        violation = 0;
    end
end






end

