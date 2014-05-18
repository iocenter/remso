function [firstOptDualApprox,errorSum] = multiplierFreeApproxs(objPartials,ax,mu,av,d,withAlgs)
% \lambda^{\top} * \frac{\partial c}{\partial \theta}
% only range space solution considered, nullspace part multiplies to 0!

firstOptDualApprox =  - sum([...
    sum(cellfun(@(Jz,dz)Jz*dz,objPartials.Jx',ax));...
    sum(cellfun(@(dz,uz,lz)((uz-lz)'*dz),ax,mu.ub.x,mu.lb.x))
    ]);

if withAlgs
    firstOptDualApprox = firstOptDualApprox - sum([...
        sum(cellfun(@(Jz,dz)Jz*dz,objPartials.Jv',av));...
        sum(cellfun(@(dz,uz,lz)((uz-lz)'*dz),av,mu.ub.v,mu.lb.v))...
        ]);
end

errorSum = sum(cellfun(@(di)sum(abs(di)),d));
end