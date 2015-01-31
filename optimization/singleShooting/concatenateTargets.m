function [ target ] = concatenateTargets(target1,target2,outDims)

% creates a function target that is the concatenation of target1 and target2


target = cell(1,numel(target1));


for k = 1:numel(target1)
    target{k} = arroba(@vertcatFunctions,[1,2,3],{target1{k},target2{k},outDims},true);
end



end

function [f,Jac]  = vertcatFunctions(xsk,uk,vsk,target1,target2,outDims,varargin)

opt = struct('partials',false,'leftSeed',[],'xRightSeeds',[],'uRightSeeds',[],'vRightSeeds',[]);
opt = merge_options(opt, varargin{:});


leftSeed1 = [];
leftSeed2 = [];

if  ~(size(opt.leftSeed,2)==0)
    leftSeed1 = opt.leftSeed(:,1:outDims(1));
    leftSeed2 = opt.leftSeed(:,outDims(1)+1:end);
else
    
end


[fK1,JacK1] = callArroba(target1,{xsk,uk,vsk},...
    'partials',opt.partials,...
    'leftSeed',leftSeed1,...
    'xRightSeeds',opt.xRightSeeds,...
    'uRightSeeds',opt.uRightSeeds,...
    'vRightSeeds',opt.vRightSeeds);

[fK2,JacK2] = callArroba(target2,{xsk,uk,vsk},...
    'partials',opt.partials,...
    'leftSeed',leftSeed2,...
    'xRightSeeds',opt.xRightSeeds,...
    'uRightSeeds',opt.uRightSeeds,...
    'vRightSeeds',opt.vRightSeeds);

f = [fK1;
    fK2];


Jac = [];
if opt.partials
    if ~(size(opt.xRightSeeds,1)==0)
        
        Jac.J = [JacK1.J;
            JacK2.J];
        
    elseif ~(size(opt.leftSeed,2)==0)
        
        Jac.Jx = JacK1.Jx + JacK2.Jx;
        Jac.Jv = JacK1.Jv + JacK2.Jv;
        Jac.Ju = JacK1.Ju + JacK2.Ju;
    else
        
        Jac.Jx = [JacK1.Jx ; JacK2.Jx];
        Jac.Jv = [JacK1.Jv ; JacK2.Jv];
        Jac.Ju = [JacK1.Ju ; JacK2.Ju];
        
    end
end





end
