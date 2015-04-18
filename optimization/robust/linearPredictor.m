function [ Ax,Av,As ] = linearPredictor(du,x,u,v,s,ss,simVars)

opt.etaRisk = 0.9;

[xs,vs,s2,xd,vd,sd,ax,Ax,av,Av,as,As]  = condensing_R(x,u,v,s,ss,'simVars',simVars,...
    'uRightSeeds',du,...
    'computeCorrection',false,...
    'computeNullSpace',true);


end

