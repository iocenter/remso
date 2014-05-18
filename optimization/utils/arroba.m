function [ arrobaFun ] = arroba(fun,argsOrder,fixedArgs,allowAdditionalArgs)
%
%  Replacement for anonymous functions. anonymous functions save too much
%  information that is not needed in our application
% 
%  See also:
%
%  callArroba

if nargin < 3
    fixedArgs  = {};
end
if nargin < 4
    allowAdditionalArgs  = false;
end

arrobaFun.f = fun;
arrobaFun.argsOrder = argsOrder;
arrobaFun.fixedArgs = fixedArgs;
arrobaFun.allowAdditionalArgs = allowAdditionalArgs;


end

