function [ s,Jac ] = outputRisks(o,varargin)

opt = struct('eta',0,'partials',false,'leftSeed',[],'oRightSeed',[]);
opt = merge_options(opt, varargin{:});

noT = numel(o{1});
eta = opt.eta;
if numel(eta) == 1 
    eta = repmat(eta,noT,1);
end

[c,jacC] = arrayfun(@(no) applyRiskMeasure(cellfun(@(or)or(no),o),eta(no)),(1:noT)','UniformOutput',false);
s = cell2mat(c);

Jac = [];
if opt.partials
    
    oDims = cellfun(@numel,o);
    pos = cumsum([0;oDims]);
    
    JacO = cellfun(@(jacCi,noi)calcOutputJacobian(jacCi,pos,noi),jacC,num2cell(1:numel(jacC))','UniformOutput',false) ;
    JacO = cell2mat(JacO);
    
    
    if  iscell(opt.oRightSeed) && (size(opt.oRightSeed{1},1)~=0)
        Jac.J = JacO*cell2mat(opt.oRightSeed);
        
    else
        oDims = cellfun(@numel,o);
        
        if  (size(opt.leftSeed,2)~=0)
            JacO = opt.leftSeed*JacO;
            Jac.Jo  = mat2cell(JacO,size(opt.leftSeed,1),oDims);
        else
            no = numel(s);
            Jac.Jo  = mat2cell(JacO,no,oDims);
        end
        Jac.Jo = Jac.Jo';
    end
end


end

function jac = calcOutputJacobian(JacCo,pos,no)
% i must be a vector of 1'
[i,j,v] = find(JacCo);

newJ = arrayfun(@(ji)pos(ji)+no,j);

jac = sparse(i,newJ,v,1,pos(end));


end

function [s,jac] = applyRiskMeasure(oi,eta)

if eta == 1  %% worst case scenario
    [s,i] = max(oi);
    jac = sparse(1,i,1,1,numel(oi));

elseif eta == 0 % mean value
    noi = numel(oi);
    [s] = mean(oi);
    jac = repmat(1/noi,1,noi);
else  % cvar with a certain value
    [s,jac] = cvar(oi,eta);
end

end
    