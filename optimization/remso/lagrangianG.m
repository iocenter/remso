function [lagG] = lagrangianG(  u,x,v,lambdaX,lambdaV,muU,muX,muV,obj,ss,varargin)

opt = struct('xs',[],'vs',[],'simVars',[]);
opt = merge_options(opt, varargin{:});

lagG = [];

[~,objPartials] = obj(x,u,v,'gradients',true);


[xs,vs,Jac,convergence,simVars,usliced] = simulateSystem(x,u,ss,'gradients',true,'guessX',opt.xs,'guessV',opt.vs,'xLeftSeed',lambdaX,'vLeftSeed',lambdaV,'simVars',opt.simVars);

gineq.Ju = muU;
gineq.Jx = muX;
gineq.Jv = muV;

geq.Jx = cellfun(@(lrs,l)lrs-l,Jac.Jx,lambdaX,'uniformOutput',false);
geq.Jv = cellfun(@(l)-l,lambdaV,'uniformOutput',false);
geq.Ju = cellfun(@(lrs)lrs,Jac.Ju,'uniformOutput',false);

lagG.Ju = cellfun(@(o,eq,ineq)o+eq+ineq,objPartials.Ju,geq.Ju,gineq.Ju,'UniformOutput',false);
lagG.Jx = cellfun(@(o,eq,ineq)o+eq+ineq,objPartials.Jx,geq.Jx,gineq.Jx,'UniformOutput',false);
lagG.Jv = cellfun(@(o,eq,ineq)o+eq+ineq,objPartials.Jv,geq.Jv,gineq.Jv,'UniformOutput',false);


end


