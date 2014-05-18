function [f,gradU,targetPartials ] = targetGrad(xs,u,vs,target,Ax,Av,ci,varargin)
% Compute the value and gradient of a target function.
%
% SYNOPSIS:
%  [f,gradU,targetPartials] = targetGrad(xs,u,vs,target,Ax,Av,ci)
%  [f,gradU,targetPartials] = targetGrad(x,u,v,ss, 'pn', pv, ...)
%
% PARAMETERS:
%
%   xs - Simulated state output.
%
%   u - cellarray containing the controls for each control
%       period.
%
%   vs - Simulated algebraic state output.
%
%   target -  General target function
%
%   Ax - gradient of the states w.r.t. the controls
%
%   Av - gradient of the algebraic states w.r.t. the controls
%
%   ci - control incidence function
%
%   'pn'/pv - List of 'key'/value pairs defining optional parameters. The
%             supported options are:
%
%   usliced - controls for each shooting interval.
%
%
% RETURNS:
%
%   f - value of the target function.
%
%   gradU - target function gradient w.r.t. u.
%
%   targetPartials - Partial derivatives of the target function
%
% SEE ALSO:
%
%


opt = struct('usliced',[]);
opt = merge_options(opt, varargin{:});

[f,targetPartials] = target(xs,u,vs,'gradients',true,'usliced',opt.usliced);

gradU = targetPartials.Ju;

for i = 1:size(Ax,1)  
    jF = ci(i);  %% lower triangular!
    for j = 1:jF
        if ~isempty(Ax{i,j})
            gradU{j} = gradU{j} + targetPartials.Jx{i} * Ax{i,j};
        end
        if ~isempty(Av{i,j})
            gradU{j} = gradU{j} + targetPartials.Jv{i} * Av{i,j};
        end
    end
end

end


