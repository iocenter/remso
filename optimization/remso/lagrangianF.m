function [ lagF,lagG] = lagrangianF(  u,x,v,lambdaX,lambdaV,muU,muX,muV,obj,ss,lbx,lbv,lbu,ubx,ubv,ubu,varargin)

opt = struct('gradients',false,'xs',[],'vs',[],'simVars',[],'withAlgs',false);
opt = merge_options(opt, varargin{:});

withAlgs = opt.withAlgs;

lagG = [];

[f,objPartials] = obj(x,u,v,'gradients',opt.gradients);


[xs,vs,Jac,convergence,simVars,usliced] = simulateSystem(x,u,ss,'gradients',opt.gradients,'guessX',opt.xs,'guessV',opt.vs,'xLeftSeed',lambdaX,'vLeftSeed',lambdaV,'simVars',opt.simVars,'withAlgs',withAlgs);

feq = sum(cellfun(@(l,xsi,xi)l*(xsi-xi),[lambdaX,lambdaV]',[xs;vs],[x;v]));

[fineq,gineq] = ineqFun(u,x,v,muU,muX,muV,lbx,lbv,lbu,ubx,ubv,ubu,'gradients',opt.gradients);

lagF = f + feq + fineq;


if opt.gradients
    
    geq.Jx = cellfun(@(lrs,l)lrs-l,Jac.Jx,lambdaX,'uniformOutput',false);
    geq.Jv = cellfun(@(l)-l,lambdaV,'uniformOutput',false);
    geq.Ju = cellfun(@(lrs)lrs,Jac.Ju,'uniformOutput',false);

    lagG.Ju = cellfun(@(o,eq,ineq)o+eq+ineq,objPartials.Ju,geq.Ju,gineq.Ju,'UniformOutput',false);
    lagG.Jx = cellfun(@(o,eq,ineq)o+eq+ineq,objPartials.Jx,geq.Jx,gineq.Jx,'UniformOutput',false);
    lagG.Jv = cellfun(@(o,eq,ineq)o+eq+ineq,objPartials.Jv,geq.Jv,gineq.Jv,'UniformOutput',false);

end



end



function [fineq,gineq] = ineqFun(u,x,v,muU,muX,muV,lbx,lbv,lbu,ubx,ubv,ubu,varargin)

opt = struct('gradients',false,'xs',[],'vs',[]);
opt = merge_options(opt, varargin{:});

gineq = [];

uCellDim = numel(u);
xCellDim = numel(x);
vCellDim = numel(v);


mu = cell2mat([muU,muX,muV]);

lb = cell2mat([lbu;lbx;lbv]);

ub = cell2mat([ubu;ubx;ubv]);

z = [u;x;v];
zDim = cellfun(@(xi)numel(xi),z);
z = cell2mat(z);


muPlus =  max(mu,0);
muMinus = min(mu,0);

fineq = muPlus*(z-ub) + muMinus*(z-lb);

if opt.gradients

    gineqMat = mu;

    gineqCell = mat2cell(gineqMat,1,zDim');

    gineq.Ju = gineqCell(1:uCellDim);
    gineq.Jx = gineqCell(uCellDim+1:uCellDim+xCellDim);
    gineq.Jv = gineqCell(uCellDim+xCellDim+1:uCellDim+xCellDim+vCellDim);

end

end