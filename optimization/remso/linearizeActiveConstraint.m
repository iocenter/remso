function [ Aact ] = linearizeActiveConstraint(activeSet,u,x,v,ss,simVars,uDims,withAlgs )


[~,JacAct ] = activeSet2TargetXV(uDims,activeSet);


[~,Aact,~,~,~,~] = simulateSystemZ(u,x,v,ss,[],'simVars',simVars,'JacTar',JacAct,'withAlgs',withAlgs);



end

