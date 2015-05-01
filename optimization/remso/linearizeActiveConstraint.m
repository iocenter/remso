function [ Aact ] = linearizeActiveConstraint(activeSet,u,x,v,ss,simVars,withAlgs )


[~,JacAct ] = activeSet2TargetXV(activeSet);


[~,Aact] = simulateSystemZ(u,x,v,ss,[],'simVars',simVars,'JacTar',JacAct,'withAlgs',withAlgs);



end

