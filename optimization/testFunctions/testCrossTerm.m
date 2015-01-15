function [ e ] = testCrossTerm( x,v,u,obj,ss,muX,muV,muU )
%{
There are 2 ways to compute the cross term by FD according to

C. Schmid, “Reduced Hessian successive programming for large-scale process optimization,” Carnegie Mellon, 1994.

Essentially, both approximates the cross term with an error of the order of
the range space solution. However, there is an important remark:  one of
the approximations is dependent on all the dual variables (
which makes a lot of sense), while the other is independent of the lagrange
multipliers of the equality constraints (which does not make too much
sense).



%}

e = [];

withAlgs = (ss.nv >0);

[f,objPartials] = obj(x,u,v,'gradients',true);
[xs,vs,Jac,convergence,simVars,usliced] = simulateSystem(x,u,ss);

[xs,vs,xd,vd,ax,Ax,av,Av] = condensing(x,u,v,ss,'computeCorrection',true,'simVars',simVars);


% gbar = g+nu
gbar.Jx =  cellfun(@(Jz,mi)(Jz+mi),objPartials.Jx,muX,'UniformOutput',false);
gbar.Ju =  cellfun(@(Jz,mi)(Jz+mi),objPartials.Ju,muU,'UniformOutput',false);
if withAlgs
    gbar.Jv = cellfun(@(Jz,mi)(Jz+mi),objPartials.Jv,muV,'UniformOutput',false);
end
[~,~,~,~,~,~,~,lambdaX,lambdaV]= simulateSystemZ(u,xd,vd,ss,[],'gradients',true,'guessX',xs,'guessV',vs,'simVars',simVars,'JacTar',gbar);


[ hhxvu ] = buildFullHessian( x,v,u,obj,ss,lambdaX,lambdaV);

uDims = cellfun(@(ui)numel(ui),u);
xDims = cellfun(@(ui)numel(ui),x);
vDims = cellfun(@(ui)numel(ui),v);



%nullspace and range space

Z = [cell2matFill(Ax);cell2matFill(Av);eye(sum(uDims))];
Y = [speye(sum(xDims+vDims));zeros(sum(uDims),sum([xDims;vDims]))];


mudx = cellfun(@(z)z',muX','UniformOutput',false);
mudv = cellfun(@(z)z',muV','UniformOutput',false);
mudu = cellfun(@(z)z',muU','UniformOutput',false);


pert = 0.00001;
lbxH = cellfun(@(zi)zi-pert,x,'UniformOutput',false);
ubxH = cellfun(@(zi)zi+pert,x,'UniformOutput',false);

lbvH = cellfun(@(zi)zi-pert,v,'UniformOutput',false);
ubvH = cellfun(@(zi)zi+pert,v,'UniformOutput',false);


% compute cross-term



if withAlgs
    gbarZ = vectorTimesZ(gbar.Jx,gbar.Ju,gbar.Jv,Ax,Av,ss.ci );
else
    gbarZ = vectorTimesZ(gbar.Jx,gbar.Ju,[],Ax,[],ss.ci );
end
[w6,stepY] = computeCrossTerm(x,u,v,ax,av,gbarZ,ss,obj,mudx,mudu,mudv,lbxH,lbvH,ubxH,ubvH,withAlgs,'xs',xs,'vs',vs);

% Cross term by definition (using FD approx for the hessian)
wF = Z'*cell2mat(hhxvu)*Y*cell2mat([ax;av]);

e = [e norm(cell2mat(w6)-wF')];


%additional check for the cross term
[ lagG] = lagrangianG( u,x,v,lambdaX,lambdaV,muU,muX,muV,obj,ss,'simVars',simVars);

xR = cellfun(@(z,dz)z+dz,x,ax,'UniformOutput',false);
vR = cellfun(@(z,dz)z+dz,v,av,'UniformOutput',false);
[xsR,vsR,~,convergenceR,simVarsR,uslicedR] = simulateSystem(xR,u,ss,'guessV',xR,'guessX',vR);
[lagGRC] = lagrangianG( u,xR,vR,lambdaX,lambdaV,muU,muX,muV,obj,ss,'simVars',simVarsR);


diffGradLag.Jx = cellfun(@(lR,l)lR-l,lagGRC.Jx,lagG.Jx,'UniformOutput',false); 
diffGradLag.Ju = cellfun(@(lR,l)lR-l,lagGRC.Ju,lagG.Ju,'UniformOutput',false);
if withAlgs
    diffGradLag.Jv = cellfun(@(lR,l)lR-l,lagGRC.Jv,lagG.Jv,'UniformOutput',false); 
end

[~,w2,convergedR,~,~,~,~ ] = simulateSystemZ(u,xd,vd,ss,[],'gradients',true,'guessV',xs,'guessX',vs,'simVars',simVars,'JacTar',diffGradLag);


e = [e norm(cell2mat(w2)-wF')];



end

