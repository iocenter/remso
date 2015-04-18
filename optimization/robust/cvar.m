function [ cVar,cVarJAC,VaR ] = cvar( c,eta )
% Computes the conditional value at risk (or average value at risk) of a
% sampled random variable c



[cVar,cVarJAC,VaR] = solveLPSort(c,eta);
%[cVar2,cVarJAC2,VaR2] = solveCPLEX(c,eta);

end

function [cVar,cVarJAC,VaR] = solveLPSort(c,eta)
% The dual problem is very easy!... lets exploit it!
% if you want to see how to solve with an LP look the function below

nr = numel(c);
nnZ = nr*(1-eta);

jz = 1/nnZ;


[~,sortIndex] = sort(c,'descend');                                                  
maxIndex = sortIndex(1:floor(nnZ));


cVarJAC = sparse(1,maxIndex,jz,1,nr,numel(maxIndex)+1);

lastV = nnZ - floor(nnZ);
if lastV > 0
	maxIndexLast = sortIndex(ceil(nnZ));
    cVarJAC(maxIndexLast) = lastV/nnZ;
else
    maxIndexLast = [];
end

activeSet = [maxIndex;maxIndexLast];


cVar = cVarJAC(activeSet)*c(activeSet);
VaR = min(c(activeSet));

end


function [cVar,cVarJAC,VaR] = solveCPLEX(c,eta)


nr = numel(c);

P = Cplex('LP - CVAR');
P.DisplayFunc = [];
P.addCols([1;repmat(1/((1-eta)*nr),nr,1)],[],[-inf;zeros(nr,1)],inf(nr+1,1));
P.addRows(c,[ones(nr,1),speye(nr)],inf(nr,1));
P.solve();

% norm(P.Model.obj' - P.Solution.dual' * P.Model.A - P.Solution.reducedcost')


cVar = full(P.Solution.objval);
cVarJAC = sparse(P.Solution.dual');
VaR = P.Solution.x(1);


end