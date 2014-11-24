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
    
    builtArgs = fun.fixedArgs;
    
    for k = 1:numel(fun.argsOrder)
        if fun.argsOrder(k)>0  %% only positive indexes are considered
            builtArgs = [builtArgs(1:(fun.argsOrder(k)-1)),...
                args{k},...
                builtArgs(fun.argsOrder(k):end)];
        end
    end
    
    
    varargout = cell(1:nargout);
    
    if isstruct(fun.f)  %% allow recursive calling of arroba functions
        [varargout{1:nargout}] = callArroba(fun.f,builtArgs,varargin{:});    
    else
        [varargout{1:nargout}] = fun.f(builtArgs{:},varargin{:});
    end
    
else  %% should be a function handle!
    varargout = cell(1:nargout);
    [varargout{1:nargout}] = fun(args{:},varargin{:});
end


end

