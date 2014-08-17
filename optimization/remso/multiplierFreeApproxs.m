function [firstOptDualApprox,errorSum] = multiplierFreeApproxs(objPartials,ax,av,xd,vd,mu,withAlgs)
% \lambda^{\top} * \frac{\partial c}{\partial \theta}
% only range space solution considered, nullspace part multiplies to 0!



xJN = cellfun(@(x,y)x-y,ax,xd,'UniformOutput',false);

firstOptDualApprox = - sum([...
    sum(cellfun(@(Jz,dz)Jz*dz,objPartials.Jx',xJN));...
    sum(cellfun(@(dz,uz,lz)((uz-lz)'*dz),xJN,mu.ub.x,mu.lb.x))
    ]);

if withAlgs
    vJN = cellfun(@(x,y)x-y,av,vd,'UniformOutput',false);
    
    firstOptDualApprox = firstOptDualApprox - sum([...
        sum(cellfun(@(Jz,dz)Jz*dz,objPartials.Jv',vJN));...
        sum(cellfun(@(dz,uz,lz)((uz-lz)'*dz),vJN,mu.ub.v,mu.lb.v))...
        ]);
end

d = [xd;vd];
errorSum = sum(cellfun(@(di)sum(abs(di)),d));


%{
% both should give same results!


duV = cell2mat(du);
xJ = cellfun(@(x,y)x-y,dx,xd,'UniformOutput',false);

firstOptDualApprox =  -duV'*M*duV  - sum([...
    sum(cellfun(@(Jz,dz)Jz*dz,objPartials.Jx',xJ));...
    sum(cellfun(@(Jz,dz)Jz*dz,objPartials.Ju',du));...
    sum(cellfun(@(dz,uz,lz)((uz-lz)'*dz),xJ,mu.ub.x,mu.lb.x))
    sum(cellfun(@(dz,uz,lz)((uz-lz)'*dz),du,mu.ub.u,mu.lb.u))...
    ]);

if withAlgs
    vJ = cellfun(@(x,y)x-y,dv,vd,'UniformOutput',false);
    
    firstOptDualApprox = firstOptDualApprox - sum([...
        sum(cellfun(@(Jz,dz)Jz*dz,objPartials.Jv',vJ));...
        sum(cellfun(@(dz,uz,lz)((uz-lz)'*dz),vJ,mu.ub.v,mu.lb.v))...
        ]);
end

d = [xd;vd];
errorSum = sum(cellfun(@(di)sum(abs(di)),d));



%}




end