function [ f,Jac ] = lastAlg( x,u,v,varargin)


opt = struct('partials',false,'leftSeed',[],'xRightSeeds',[],'uRightSeeds',[],'vRightSeeds',[]);
opt = merge_options(opt, varargin{:});

f = v(end);

Jac = [];
if opt.partials
    if ~isempty(opt.vRightSeeds)
        Jac.J = opt.vRightSeeds(end,:);
    else
        nv = numel(v);
        Jv = sparse(1,nv,1,1,nv);
        if size(opt.leftSeed,2)>0
            Jac.Ju = sparse(size(opt.leftSeed,1),numel(u));
            Jac.Jx = sparse(size(opt.leftSeed,1),numel(x));
            Jac.Jv = opt.leftSeed*Jv;
        else
            Jac.Ju = sparse(1,numel(u));
            Jac.Jx = sparse(1,numel(x));
            Jac.Jv = Jv;
        end
        
    end
end
end

