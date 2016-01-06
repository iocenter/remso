function [ lowActive,upActive ] = checkConstraintActivity(lb,ub,var,varargin)
%  lowActive = (ldx + tol >= dx)
%  lowActive = (udx - tol <= dx)


opt = struct('tol',1e-5);
opt = merge_options(opt, varargin{:});

upActive  = cellfun(@(upBound,var)  upBound-opt.tol < var,ub,var,'UniformOutput',false);
lowActive = cellfun(@(lowBound,var)lowBound+opt.tol > var,lb,var,'UniformOutput',false);

end

