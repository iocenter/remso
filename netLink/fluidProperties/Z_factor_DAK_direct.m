function Z = Z_factor_DAK_direct(p,rho_g_sc,T_abs)
% Z = Z_factor_DAK_direct(p,rho_g_sc,T_abs)
%
% Calculates the Z-factor for a given reduced pressure and reduced temperature.
% Use is made of the correlation of Dranchuk & Abu-Kasem (1975) to approximate the 
% Standing & Katz (1942) chart. Equal to Z_factor_DAK.m but takes different arguments to
% avoid the need to compute intermediate steps (pseudo properties). 
%
% The range of validity for the approximation is
% 0.2 < p_pr < 30 and 1.0 < T_pr < 3.0 .
%
% Z = Z factor, -
% p = pressure, Pa
% rho_g_sc = gas density at standard conditions, kg/m.^3
% T_abs = absolute temperature, K
%
% JDJ, 27-10-13, last revised 27-10-13
%%TODO: include initial guess for zfactor

    % Compute pseudo-reduced properties:
    p_pc = pres_pseu_crit_Sutton(rho_g_sc);
    T_pc = temp_pseu_crit_Sutton(rho_g_sc);
    p_pr = p ./ p_pc;
    T_pr = T_abs ./ T_pc;

    % Check validity of input:
    % Initiate Z with the
    if p_pr > 30
        p_pr 
        warning('Warning: Pseudo-reduced pressure too large (i.e. above 30) in Z_factor_DAK.')
    end
    if T_pr > 3.0
        T_pr
        warning('Warning: Pseudo-reduced temperature too large (i.e. above 3.0) in Z_factor_DAK.')
    end

    % Coefficients:
    a1 =  0.3265;
    a2 = -1.0700;
    a3 = -0.5339;
    a4 =  0.01569;
    a5 = -0.05165;
    a6 =  0.5475;
    a7 = -0.7361;
    a8 =  0.1844;
    a9 =  0.1056;
    a10 = 0.6134;
    a11 = 0.7210;

    c = 0.27 * p_pr./T_pr;

    b1 = c .* (a1 + a2./T_pr + a3./T_pr.^3 + a4./T_pr.^4 + a5./T_pr.^5);
    b2 = c.^2 .* (a6 + a7./T_pr + a8./T_pr.^2);
    b3 = c.^5 .* a9.*(a7./T_pr + a8./T_pr.^2);
    b4 = c.^2 .* a10./T_pr.^3;
    b5 = c.^2 .* a11;
    b6 = b4 .* b5;
    
    %% create doubles for ADI parameters
    double_b1 = double(b1);
    double_b2 = double(b2);
    double_b3 = double(b3);
    double_b4 = double(b4);
    double_b5 = double(b5);
    double_b6 = double(b6);

    % Initiate Z with the Papay correlation:
    Z_0 = Z_factor_Papay(double(p_pr),double(T_pr));
    Z = Z_0;

    % Improve the result with Newton Raphson iteration:
    tol_abs = 1.e-6; % Absolute convergence criterion
    tol_rel = 1.e-9; % Relative convergence criterion
    max_iter = 1000; % Maximum allowed number of iterations 
    max_diff = 0.3; % Maximum allowed absolute difference in Z per iteration step
    iter = 0; % Iteration counter
    repeat = 1; 
    
    cond_convergence = true(numel(double(p_pr)),1);
    absDiff = zeros.*cond_convergence;
    rel_diff = zeros.*cond_convergence;
    Z_old = Z;
    while repeat > 0
       if iter > max_iter
          p_pr
          T_pr
          Z_0
          Z
          error('Error: Maximum allowed number of iterations exceeded in Z_factor_DAK.')
       end   
       iter = iter+1;
       Z_old(cond_convergence) = Z(cond_convergence);

       [fZ, dfZdZ] = fzCalc(Z_old(cond_convergence), double_b1(cond_convergence), double_b2(cond_convergence), double_b3(cond_convergence), double_b4(cond_convergence), double_b5(cond_convergence), double_b6(cond_convergence), true);    

       Z(cond_convergence) = Z_old(cond_convergence) - fZ(cond_convergence)./dfZdZ(cond_convergence); % Newton Raphson iteration
       absDiff(cond_convergence) = Z(cond_convergence)-Z_old(cond_convergence);
     
       
       cond_maxdiff = abs(absDiff(cond_convergence)) > max_diff;
       if any(cond_maxdiff) % Check if steps are too large
          Z(cond_convergence & cond_maxdiff) = Z_old(cond_convergence & cond_maxdiff) + max_diff* sign(absDiff(cond_convergence & cond_maxdiff)); % Newton Raphson iteration with reduced step size
          absDiff(cond_convergence & cond_maxdiff) = max_diff;   
       end       
       
       rel_diff(cond_convergence) = absDiff(cond_convergence)./Z_old(cond_convergence);        
       
       cond_rel = abs(rel_diff(cond_convergence)) > tol_rel;      
       cond_abs = abs(absDiff(cond_convergence)) > tol_abs;     
       
       if any(cond_abs) % Check for convergence
          repeat = 1;
       else
          if any(cond_rel)
             repeat = 1;
          else
             repeat = 0;
          end
       end   
       cond_convergence = cond_abs | cond_rel;       
    end     
    
    if isa(p, 'ADI')
        [fz, dfzdz] = fzCalc(Z_old, b1, b2, b3, b4, b5, b6, true);
        for jp = 1:numel(fz.jac)
            fz.jac{jp} = - diag(1./dfzdz)*fz.jac{jp};
        end
        Z = ADI(Z,fz.jac);
        
    end
    
end

function [fz, dfzdz] = fzCalc(Z_old, b1, b2, b3, b4, b5, b6, gradFlag)    
    if nargin < 8
        gradFlag = true;
    end;
    
    help01 = Z_old - b1.*Z_old.^-1 - b2.*Z_old.^-2 + b3.*Z_old.^-5;
    help02 = -(b4.*Z_old.^-2 + b6.*Z_old.^-4) .* exp(-b5.*Z_old.^-2) - 1;
    fz = help01 + help02;

    if gradFlag
         help03 = 1 + double(b1).*Z_old.^-2 + 2.*double(b2).*Z_old.^-3 - 5.*double(b3).*Z_old.^-6;
         help04 = (2.*double(b4).*Z_old.^-3 - 2.*double(b4).*double(b5).*Z_old.^-5 + 4.*double(b6).*Z_old.^-5 - 2.*double(b5).*double(b6).*Z_old.^-7)...
        .* exp(-double(b5).*Z_old.^-2);
        dfzdz = help03 + help04;
    else
        dfzdz = [];
    end
    
end


