function [ zfac ] = gasZFactor(yg, t, p)
%      t  - oC
%      p  - barsa
%      yg - spesific gravity

t = convtemp(t,'K','C');
p = p./barsa;

%   Calculating pseudocritical properties
t_pc = 169.2 +349.5.*yg - 74.*yg.^2;
p_pc = 756.8 - 131.*yg - 3.6.*yg.^2;


%    calculating pseudo reduced properties
t_pr = (t + 460)./t_pc;
p_pr = p./p_pc;

t = 1./t_pr;
a = 0.06125.*t.*exp(-1.2.*(1-t).^2);

y = 0.001;

% Improve the result with Newton Raphson iteration:
i = 0;
tol_abs = 1.e-6; % Absolute convergence criterion
tol_rel = 1.e-9; % Relative convergence criterion
max_iter = 1000; % Maximum allowed number of iterations
max_diff = 0.3; % Maximum allowed absolute difference in Z per iteration step
repeat = 1;

fy = 1;

while  (repeat && (i < max_iter))
    y_old = y;
    
    fy = -a .* p_pr + (y + y.^2 + y.^3 - y.^4)./ (1-y).^3 - (14.76.* t - 9.76.* t.^2 + 4.58.*t.^3) .* y.^2 + (90.7.*t - 242.2.*t.^2 + 42.4 .*t.^3).*y.^(2.18 + 2.82.*t);
    
    dfY = (1 + 4.*y + 4.*y.^2 - 4.*y.^3 + y.^4) ./ (1-y).^4 - (29.52 .* t - 19.52 .* t.^2 + 9.16 .* t.^3) .* y + (2.18 + 2.82 .* t) .* (90.7 .* t - 242.2 .* t.^2 + 42.4 .* t.^3) .* y.^(1.18 + 2.82 .* t);
    
    y = y - fy ./ dfY;
    
    
    diff = y-y_old;
    if abs(diff) > max_diff % Check if steps are too large
        y = y_old + max_diff * sign(diff); % Newton Raphson iteration with reduced step size
        diff = max_diff;
    end
    rel_diff = diff./y_old;
    if abs(diff) > tol_abs % Check for convergence
        repeat = 1;
    else
        if abs(rel_diff) > tol_rel
            repeat = 1;
        else
            repeat = 0;
        end
    end
    i = i+1;
end

zfac = a.*p_pr./y;

end

% function [fy, dfy] = estimateZfactor(t, y, a, p_pr)
%     fy = -a .* p_pr + (y + y.^2 + y.^3 - y.^4)./(1 - y).^3 - (14.76 .* t - 9.76 .* t.^2 + 4.58 .* t.^3) .* y.^2 + (90.7 .* t - 242.2 .* t.^2 + 42.4 .* t.^3) .*y.^(2.18 + 2.82 .* t);
% 
%     dfy = (1 + 4 .* y + 4 .* y.^2   - 4 .* y.^3   + y.^4)./(1 - y).^4 - (29.52 .* t - 19.52 .*   t.^2   + 9.16 .*t.^3)  .* y + (2.18 + 2.82 .* t) .* (90.7 .* t - 242.2 .*   t.^2  + 42.4 .*   t.^3 ) .*   y.^(1.18 + 2.82 .* t);
% end


