function f = Moody_friction_factor(epsilon,N_Re)
% function f = Moody_friction_factor(epsilon,N_Re)
%
% Computes the friction factor for pipe flow according to the Moody (1944) diagram. In the
% turbulent region, the implicit Colebrook (1939) expression is used to compute the friction
% factor iteratively via subsequent substitution.
%
% epsilon = dimensionless roughness, -
% f = friction factor, -
% N_Re = Reynolds number, -
%
% JDJ, 12-02-02, last revised 10-05-13
%{
Changes by Thiago and Codas
Make the function compatible with ADI objects
Gradient correction regarding the iterative procedure
%}
if N_Re < 2000 % laminar regime
    f = 64./N_Re;
else % turbulent or transitional regime
    if N_Re < 3000 % transitional regime: prepare for interpolation
        f_lam_max = 64./2000; % highest laminar value
        alpha = (N_Re-2000)./(3000-2000); % interpolation parameter
        N_Re_work = 3000; % sets N_Re_work to compute lowest turbulent value
    else % turbulent regime
        N_Re_work = N_Re;
    end
    % Initialize f_work with the Zigrang and Sylvester (1985) approximation for the Colebrook
    % (1939) friction factor:
    help01 = 2*epsilon./3.7 + 13./N_Re_work;
    help02 = (5.02./N_Re_work)*log(help01)/log(10);
    f_work = 1./(-2*log(2*epsilon./3.7 - help02)./log(10)).^2;

    % Improve the result through iteration:
    tol_abs = 1.e-9; % absolute convergence criterion
    tol_rel = 1.e-8; % relative convergence criterion
    max_iter = 100; % maximum allowed number of iterations
    iter = 0; % iteration counter
    repeat = 1;
    while repeat > 0
        if iter > max_iter
            error('Error: Maximum allowed number of iterations exceeded in Moody_friction_factor.')
        end
        iter = iter+1;
        f_old = double(f_work);

        % Improve the estimate:
        [f_work] = ColebrookIter(double(N_Re_work),double(f_old),double(epsilon));

        % Check for convergence:
        diff = f_work-f_old;
        rel_diff = diff./f_old;
        if abs(diff) > tol_abs
            repeat = 1;
        else
            if abs(rel_diff) > tol_rel
                repeat = 1;
            else
                repeat = 0;
            end
        end
    end
    if isa(N_Re_work,'ADI')
        assert(~isa(epsilon,'ADI'));
        
        N_Re_work2 = ADI(N_Re_work.val,[N_Re_work.jac,0]);
        f_work2 = ADI(f_work,[cellfun(@(x)x*0,N_Re_work.jac,'UniformOutput',false),1]);
        
        [f_work2] = ColebrookIter(N_Re_work2,f_work2,epsilon);
        
        
        invDfzdz = 1./(1-f_work2.jac{end});
        J = N_Re_work.jac;
        for jp = 1:numel(J)
            J{jp} = bsxfun(@times, invDfzdz, f_work2.jac{jp});
        end
        f_work = ADI(f_work,J);    
    end
    
    
    
    
    
    if N_Re < 3000 % transitional regime: interpolate between highest laminar and lowest
        % turbulent values
        f_turb_min = f_work;
        f = f_lam_max + alpha * (f_turb_min - f_lam_max);
    else % turbulent regime
        f = f_work;
    end
end

end

function [f_work] = ColebrookIter(N_Re_work,f_old,epsilon)
    
help03 = 18.7./(N_Re_work.*(f_old.^(0.5)));
f_work = 1./(1.74 - 2*(log(2*epsilon + help03)/log(10))).^2;

end


