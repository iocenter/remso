function [ targetsC,lbt,ubt,sparsity] = outputVarsBoundSelector(lbx,ubx,lbv,ubv,uDim,ci)
%
% Create a function that returns the state constraints x and v with finite
% bounds
%

K = size(lbx,1);

targetsC = cell(K,1);
lbt = cell(K,1);
ubt = cell(K,1);
sparsity = cell(K,1);

uT = callArroba(ci,{K});
uDimCum = cumsum(uDim);

for k = 1:K
    
    finiteBoundsX = or(isfinite(lbx{k}),isfinite(ubx{k}));
    finiteBoundsV = or(isfinite(lbv{k}),isfinite(ubv{k}));
    
    lbxC = lbx{k}(finiteBoundsX);
    ubxC = ubx{k}(finiteBoundsX);
    
    lbvC = lbv{k}(finiteBoundsV);
    ubvC = ubv{k}(finiteBoundsV);
    
    
    lbt{k} = [lbxC;
              lbvC];
    ubt{k} = [ubxC;
              ubvC];
    
   sparsity{k} = zeros(size(lbxC,1)+size(lbvC,1),uDimCum(uT));
   sparsity{k}(:,1:uDimCum(callArroba(ci,{k}))) = 1;       
          
          
    
    targetsC{k} = arroba(@outputSelector,[1,2,3],...
        {finiteBoundsX,...
        finiteBoundsV},...
        true);

end

end

function [T,Jac] = outputSelector(xk,uk,vk,ix,iv,varargin)

opt = struct('partials',false,'leftSeed',[],'xScale',[],'vScale',[],'uScale',[],'xRightSeeds',[],'uRightSeeds',[],'vRightSeeds',[]);
opt = merge_options(opt, varargin{:});


T = [xk(ix);
     vk(iv)];

 
Jac = [];

if opt.partials

    
    if (size(opt.leftSeed,2)==0)
        opt.leftSeed = 1;
    end
    
    nxOut = sum(ix);
    nvOut = sum(iv);
	tSize = nxOut+nvOut;
    
    xSize = numel(xk);
    vSize = numel(vk);
    uSize = numel(uk);
    
	xJx = eye(xSize); xJx = xJx(ix,:);    
    Jx = zeros(tSize,xSize);
    Jx(1:nxOut,:) = xJx;
   
    vJv = eye(vSize); vJv = vJv(iv,:);
    Jv = zeros(tSize,vSize);
    Jv(nxOut+1:end,:) = vJv;
    
	Ju = zeros(tSize,uSize);
           
    
    if (size(opt.xRightSeeds,1)==0)
        Jac.Jx = opt.leftSeed * Jx;        
        Jac.Jv = opt.leftSeed * Jv;
        Jac.Ju = opt.leftSeed * Ju;
    else
        Jac.J = opt.leftSeed*(Jx*opt.xRightSeeds + Jv*opt.vRightSeeds);  %+ Ju*opt.uRightSeeds is zero! 
    end
    
    
    
end
 
 


end