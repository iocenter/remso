function [ f,Jac ] = concatenateTargetK(k,xsk,uk,vsk,target,targetSizes,varargin)
%
%  evaluate the functions in target and concatenate the results

opt = struct('partials',false,'leftSeed',[],'xScale',[],'vScale',[],'uScale',[],'xRightSeeds',[],'uRightSeeds',[],'vRightSeeds',[]);
opt = merge_options(opt, varargin{:});


tSize = sum(targetSizes);
i1 = 1+sum(targetSizes(1:k-1));
iend = i1+targetSizes(k)-1;


if ~(size(opt.leftSeed,2)==0)
   leftSeed = opt.leftSeed(:,i1:iend);
else
   leftSeed = []; 
end


[fK,JacK] = callArroba(target,{xsk,uk,vsk},...
                    'partials',opt.partials,...
                    'leftSeed',leftSeed,...
                    'xScale',opt.xScale,...
                    'vScale',opt.vScale,...
                    'uScale',opt.uScale,...
                    'xRightSeeds',opt.xRightSeeds,...
                    'uRightSeeds',opt.uRightSeeds,...
                    'vRightSeeds',opt.vRightSeeds);
                
                
 f = zeros(tSize,1);
 f(i1:iend) = fK;
 

Jac = [];
if opt.partials
    if ~(size(opt.xRightSeeds,1)==0)
        
        Jac.J = zeros(tSize,size(opt.xRightSeeds,2));
        Jac.J(i1,iend,:) = JacK.J;
        
    else
        
        Jac.Jx = zeros(tSize,size(xsk,1));
        Jac.Jv = zeros(tSize,size(vsk,1));
        Jac.Ju = zeros(tSize,size(uk,1));

        Jac.Jx(i1:iend,:) = JacK.Jx;
        Jac.Jv(i1:iend,:) = JacK.Jv;
        Jac.Ju(i1:iend,:) = JacK.Ju;
        
    end        
end





end

