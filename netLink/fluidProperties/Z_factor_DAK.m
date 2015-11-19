function Z = Z_factor_DAK(p_pr,T_pr)
% Z = Z_factor_DAK(p_pr,T_pr)
%
% Calculates the Z-factor for a given reduced pressure and reduced temperature.
% Use is made of the correlation of Dranchuk & Abu-Kasem (1975) to approximate the 
% Standing & Katz (1942) chart.
%
% The range of validity for the approximation is
% 0.2 < p_pr < 30 and 1.0 < T_pr < 3.0 .
%
% Z = Z factor, -
% p_pr = pseudo-reduced pressure, -
% T_pr = pseudo-reduced temperature, -
%
% JDJ, 05-03-01, last revised 01-05-12

% Check validity of input:
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

c = 0.27 * p_pr/T_pr;

b1 = c * (a1 + a2/T_pr + a3/T_pr^3 + a4/T_pr^4 + a5/T_pr^5);
b2 = c^2 * (a6 + a7/T_pr + a8/T_pr^2);
b3 = c^5 * a9*(a7/T_pr + a8/T_pr^2);
b4 = c^2 * a10/T_pr^3;
b5 = c^2 * a11;
b6 = b4 * b5;

% Initiate Z with the Papay correlation:
Z_0 = Z_factor_Papay(p_pr,T_pr);
Z = Z_0;

% Improve the result with Newton Raphson iteration:
tol_abs = 1.e-6; % Absolute convergence criterion
tol_rel = 1.e-9; % Relative convergence criterion
max_iter = 1000; % Maximum allowed number of iterations 
max_diff = 0.3; % Maximum allowed absolute difference in Z per iteration step
iter = 0; % Iteration counter
repeat = 1; 
while repeat > 0
   if iter > max_iter
      p_pr
      T_pr
      Z_0
      Z
      error('Error: Maximum allowed number of iterations exceeded in Z_factor_DAK.')
   end   
   iter = iter+1;
   Z_old = Z;
   
   help01 = Z_old - b1*Z_old^-1 - b2*Z_old^-2 + b3*Z_old^-5;
   help02 = -(b4*Z_old^-2 + b6*Z_old^-4) * exp(-b5*Z_old^-2) - 1;
   fZ = help01 + help02;
   
   help03 = 1 + b1*Z_old^-2 + 2*b2*Z_old^-3 - 5*b3*Z_old^-6;
   help04 = (2*b4*Z_old^-3 - 2*b4*b5*Z_old^-5 + 4*b6*Z_old^-5 - 2*b5*b6*Z_old^-7)...
       * exp(-b5*Z_old^-2);
   dfZdZ = help03 + help04;
   
   Z = Z_old - fZ/dfZdZ; % Newton Raphson iteration
   diff = Z-Z_old;
   if abs(diff) > max_diff % Check if steps are too large
      Z = Z_old + max_diff * sign(diff); % Newton Raphson iteration with reduced step size
      diff = max_diff;   
   end
   rel_diff = diff/Z_old;
   if abs(diff) > tol_abs % Check for convergence
      repeat = 1;
   else
      if abs(rel_diff) > tol_rel
         repeat = 1;
      else
         repeat = 0;
      end
   end
end 