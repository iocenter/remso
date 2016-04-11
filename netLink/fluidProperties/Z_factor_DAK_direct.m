function Z = Z_factor_DAK_direct(p,rho_g_sc,T_abs,zAns)
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
%{
Changes by Codas and Thiago:
Extend compatibility for ADI objects

Vectorize inputs

%%TODO: include initial guess for zfactor
%}

if nargin < 4
    zAns = [];
end


    % Compute pseudo-reduced properties:
    p_pc = pres_pseu_crit_Sutton(rho_g_sc);
    T_pc = temp_pseu_crit_Sutton(rho_g_sc);
    p_pr = p ./ p_pc;
    T_pr = T_abs ./ T_pc;

    % Check validity of input:
    % Initiate Z with the
    if any(p_pr > 30)
        p_pr 
        warning('Warning: Pseudo-reduced pressure too large (i.e. above 30) in Z_factor_DAK.')
    end
    if any(T_pr > 3.0)
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

    tol_abs = 1.e-6; % Absolute convergence criterion
    if nargin < 4
    % Initiate Z with the Papay correlation:
    Z_0 = Z_factor_Papay(double(p_pr),double(T_pr));
    Z = Z_0;

    % Improve the result with Newton Raphson iteration:
    tol_rel = 1.e-9; % Relative convergence criterion
    max_iter = 1000; % Maximum allowed number of iterations 
    max_diff = 0.3; % Maximum allowed absolute difference in Z per iteration step
    iter = 0; % Iteration counter
    repeat = 1; 
    
    cv = true(numel(double(p_pr)),1);
    absDiff = zeros(size(cv));
    rel_diff = zeros(size(cv));
    Z_old = Z;
    fZ = zeros(size(cv));
    dfZdZ = zeros(size(cv));
    while repeat > 0
       if iter > max_iter
          p_pr
          T_pr
          Z_0
          Z
          error('Error: Maximum allowed number of iterations exceeded in Z_factor_DAK.')
       end   
       iter = iter+1;
       Z_old(cv) = Z(cv);

       [fZ(cv), dfZdZ(cv)] = fzCalc(Z_old(cv), double_b1(cv), double_b2(cv), double_b3(cv), double_b4(cv), double_b5(cv), double_b6(cv), true);    

       Z(cv) = Z_old(cv) - fZ(cv)./dfZdZ(cv); % Newton Raphson iteration
       absDiff(cv) = Z(cv)-Z_old(cv);
     
       
       cond_maxdiff = abs(absDiff(cv)) > max_diff;
       if any(cond_maxdiff) % Check if steps are too large
          Z(cv & cond_maxdiff) = Z_old(cv & cond_maxdiff) + max_diff* sign(absDiff(cv & cond_maxdiff)); % Newton Raphson iteration with reduced step size
          absDiff(cv & cond_maxdiff) = max_diff;   
       end       
       
       rel_diff(cv) = absDiff(cv)./Z_old(cv);        
       
       cond_rel = abs(rel_diff(cv)) > tol_rel;      
       cond_abs = abs(absDiff(cv)) > tol_abs;     
       
       cv(cv) = cond_abs | cond_rel; 
       if any(cond_abs) % Check for convergence
          repeat = 1;
       else
          if any(cond_rel)
             repeat = 1;
          else
             repeat = 0;
             zAns = Z;
          end
       end   

    end     
    end
       [fz, dfzdz] = fzCalc(zAns, b1, b2, b3, b4, b5, b6, true);
        Z =  zAns - fz./dfzdz;  %% first order gradient correction!  Equivalent to a single newton iteration
    
        if any(abs(Z-zAns) > tol_abs)
            warning('Z correction not within tolerances');
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


