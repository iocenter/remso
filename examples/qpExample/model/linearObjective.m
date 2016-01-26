function [obj,Jac] = linearObjective(x,u,v,cx,cu,cv,varargin)

opt = struct('gradients',false,'leftSeed',[],'xRightSeeds',[],'uRightSeeds',[],'vRightSeeds',[],'usliced',[]);
opt = merge_options(opt, varargin{:});

obj = sum(cellfun(@(zi)cx*zi,x)) + sum(cellfun(@(zi)cu*zi,u)) + sum(cellfun(@(zi)cv*zi,v));

Jac = [];
if opt.gradients
    
    Jx = repmat(cx,1,numel(x));
    Ju = repmat(cu,1,numel(u));
    Jv = repmat(cv,1,numel(v));
    
    if ~(size(opt.xRightSeeds,1)==0)
        
        Jac.J = Jx*cell2mat(opt.xRightSeeds) + Ju*cell2mat(opt.uRightSeeds) + Jv*cell2mat(opt.vRightSeeds) ;
        
    else
        
        xDims = cellfun(@numel,x);
        uDims = cellfun(@numel,u);
        vDims = cellfun(@numel,v);
        
        if ~(size(opt.leftSeed,2)==0)
            
            Jac.Jx = mat2cell(opt.leftSeed*Jx,size(opt.leftSeed,1),xDims);
            Jac.Ju = mat2cell(opt.leftSeed*Ju,size(opt.leftSeed,1),uDims);
            Jac.Jv = mat2cell(opt.leftSeed*Jv,size(opt.leftSeed,1),vDims);
            
        else
            
            Jac.Jx = mat2cell(Jx,1,xDims);
            Jac.Ju = mat2cell(Ju,1,uDims);
            Jac.Jv = mat2cell(Jv,1,vDims);
        end
    end
    
end



end

