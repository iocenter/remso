function [vals] = extractCompressIneq(vals,activeVars,zDims)
%  
%  The QP problem considered by prsqpStep does not consider all
%  constraints. Therefore it is necessary to map the considered constraints
%  slacks and dual variables to the corresponding vector considering all
%  variables
%
%
%  This is helper function for this purpose.
%
%  vals - Values obtained during the QP process
%  activeVars - boolean vector of the dimension of the variables, with true
%               in all position considered within vals
%  nv - number of variables per cell block in vals.
%
activeVarsVector = cell2mat(activeVars);
valsM = zeros(size(activeVarsVector));


valsM(activeVarsVector) = vals;

vals = mat2cell(valsM,zDims,1);


end