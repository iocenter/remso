function [ varargout ] = callArroba(fun,args,varargin)
%
%  Calls a function defined by arroba or @
%
%  See also:
%
%  arroba
%


if isstruct(fun)
    assert( ((nargin <= 2) || fun.allowAdditionalArgs) )
    nDefaultArgs = numel(fun.argsOrder);
    assert( (numel(args) <= nDefaultArgs) | (fun.allowAdditionalArgs)  ) %% ups... you gave more arguments than needed!
    
    nBuiltArgs = sum(fun.argsOrder>0) + numel(fun.fixedArgs);
    builtArgs = cell(1,nBuiltArgs);
    remaingArgsIndex = 1:nBuiltArgs;    
    
    for k = 1:nDefaultArgs
        if fun.argsOrder(k)>0  %% only positive indexes are considered
            builtArgs{fun.argsOrder(k)} = args{k};
            remaingArgsIndex = setdiff(remaingArgsIndex,fun.argsOrder(k));
        end
    end
    for k = 1:numel(remaingArgsIndex)
        builtArgs{remaingArgsIndex(k)} = fun.fixedArgs{k};
    end
    
    
    varargout = cell(1:nargout);
    
    remainingArgs = args(nDefaultArgs+1:end);
    
    if isstruct(fun.f)  %% allow recursive calling of arroba functions
        [varargout{1:nargout}] = callArroba(fun.f,builtArgs,remainingArgs{:},varargin{:});    
    else
        [varargout{1:nargout}] = fun.f(builtArgs{:},remainingArgs{:},varargin{:});
    end
    
else  %% should be a function handle!
    varargout = cell(1:nargout);
    [varargout{1:nargout}] = fun(args{:},varargin{:});
end


end

