function [ s,Jac ] = outputRisks(o,varargin)

opt = struct('eta',0.9,'partials',false,'leftSeed',[],'oRightSeed',[]);
opt = merge_options(opt, varargin{:});


[c,jacC] = arrayfun(@(no) cvar(cellfun(@(si)si(no),o),opt.eta),(1:numel(o{1}))','UniformOutput',false);
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
    end
end


end

function jac = calcOutputJacobian(JacCo,pos,no)
% i must be a vector of 1'
[i,j,v] = find(JacCo);

newJ = arrayfun(@(ji)pos(ji)+no,j);

jac = sparse(i,newJ,v,1,pos(end));


end