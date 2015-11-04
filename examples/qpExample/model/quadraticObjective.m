function [obj,Jac] = quadraticObjective(x,u,v,cx,cu,cv,Qx,Qu,Qv,varargin)

opt = struct('gradients',false,'leftSeed',[],'xRightSeeds',[],'uRightSeeds',[],'vRightSeeds',[],'usliced',[],'bias',0,'scale',1);
opt = merge_options(opt, varargin{:});

obj = sum(cellfun(@(zi)cx*zi,x)) + sum(cellfun(@(zi)cu*zi,u)) + sum(cellfun(@(zi)cv*zi,v));
obj = obj + 0.5*(...
    sum(cellfun(@(zi)zi'*Qx*zi,x)) + ...
    sum(cellfun(@(zi)zi'*Qu*zi,u)) + ...
    sum(cellfun(@(zi)zi'*Qv*zi,v)));
obj = obj + opt.bias;
obj = obj * opt.scale;

Jac = [];
if opt.gradients
    
    Jx = opt.scale*(repmat(cx,1,numel(x)) + cell2mat(cellfun(@(zi)zi'*Qx,x','UniformOutput',false)));
    Ju = opt.scale*(repmat(cu,1,numel(u)) + cell2mat(cellfun(@(zi)zi'*Qu,u','UniformOutput',false)));
    Jv = opt.scale*(repmat(cv,1,numel(v)) + cell2mat(cellfun(@(zi)zi'*Qv,v','UniformOutput',false)));
    
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

