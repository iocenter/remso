function [ Aact ] = linearizeActiveConstraint(activeSet,u,x,v,ss,simVars,uDims,withAlgs )


[~,JacAct ] = activeSet2TargetXV(uDims,withAlgs,activeSet);


[~,Aact,~,~,~,~] = simulateSystemZ(u,x,v,ss,[],'simVars',simVars,'JacTar',JacAct,'withAlgs',withAlgs);



end

