function [feasible,lowActive,upActive,violation ] = checkConstraintFeasibility(dx,ldx,udx,varargin )
%
%  lowActive = (ldx + tol >= dx)
%  upActive = (udx - tol <= dx)
%
%  violation = sum of all the constraints violation (norm1 of the violation)
%
%  feasible = (violation <= opt.primalFeasTol);


opt = struct('primalFeasTol',1e-6);
opt = merge_options(opt, varargin{:});


[ lowActive,upActive ] = checkConstraintActivity(ldx,udx,dx,'tol',opt.primalFeasTol);
[ violation ] = checkConstraintViolation(dx,ldx,lowActive,udx,upActive );

feasible = (violation <= opt.primalFeasTol);


end





