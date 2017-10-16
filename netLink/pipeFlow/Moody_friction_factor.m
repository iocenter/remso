function [f,f_work] = Moody_friction_factor(epsilon,N_Re,f_work_Ans)
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
%Changes by Thiago and Codas
%Make the function compatible with ADI objects
%Gradient correction regarding the iterative procedure
%}
f_work = N_Re;    %% initialize
f = N_Re;         %% initialize
N_Re_work = N_Re; %% initialize;

c1 = N_Re < 2000;
if any(c1) % laminar regime
    f(c1) = 64./N_Re(c1);
end
if any(~c1) % turbulent or transitional regime
    c2 = N_Re < 3000 ;
    if any(~c1 & c2) % transitional regime: prepare for interpolation
        f_lam_max = 64./2000; % highest laminar value
        alpha = (N_Re(~c1 & c2)-2000)./(3000-2000); % interpolation parameter
        N_Re_work(~c1 & c2) = 3000; % sets N_Re_work to compute lowest turbulent value
    end
    if any(~c1 & ~c2) % turbulent regime
        N_Re_work(~c1 & ~c2) = N_Re(~c1 & ~c2);
    end
    % Initialize f_work with the Zigrang and Sylvester (1985) approximation for the Colebrook
    % (1939) friction factor:
    help01 = 2*epsilon(~c1)./3.7 + 13./N_Re_work(~c1);
    help02 = (5.02./N_Re_work(~c1)).*log(help01)/log(10);
    
    tol_abs = 1.e-7; % absolute convergence criterion
    if nargin <3 || isempty(f_work_Ans) || any(isnan(f_work_Ans))
        f_work(~c1) = 1./(-2*log(2*epsilon(~c1)./3.7 - help02)./log(10)).^2;
        
        % Improve the result through iteration:
        tol_rel = 1.e-8; % relative convergence criterion
        max_iter = 100; % maximum allowed number of iterations
        iter = 0; % iteration counter
        repeat = 1;
        cv = ~c1;
        f_old = double(f_work);
        f_work_Ans = f_old;
        while repeat > 0
            if iter > max_iter
                error('Error: Maximum allowed number of iterations exceeded in Moody_friction_factor.')
            end
            iter = iter+1;
            f_old(cv) = double(f_work(cv));
            
            % Improve the estimate:
            [f_work(cv)] = ColebrookIter(double(N_Re_work(cv)),double(f_old(cv)),double(epsilon(cv)),false);
            
            % Check for convergence:
            diff = f_work(cv)-f_old(cv);
            rel_diff = diff./f_old(cv);
            
            absConv = abs(diff) > tol_abs;
            relConv = abs(rel_diff) > tol_rel;
            cv(cv) = absConv | relConv;
            if any(absConv)
                repeat = 1;
            else
                if any(relConv)
                    repeat = 1;
                else
                    repeat = 0;
                    f_work_Ans(~c1) = f_work(~c1);
                end
            end
        end       
    end
    
    [f_work(~c1),dCdf]= ColebrookIter(N_Re_work(~c1),f_work_Ans(~c1),epsilon(~c1),true);
    residual = f_work(~c1)-f_work_Ans(~c1);
    dx  = - residual./(dCdf-1);
    f_work(~c1) = f_work_Ans(~c1) + dx; %% first order gradient correction!  Equivalent to a single newton iteration

    if any(abs(double(dx)) > tol_abs)
        warning('friction factor not within tolerances');
    end
        
    
    
    
    
    if any(~c1 & c2) % transitional regime: interpolate between highest laminar and lowest
        % turbulent values
        f_turb_min = f_work(~c1 & c2);
        f(~c1 & c2) = f_lam_max + alpha .* (f_turb_min - f_lam_max);
    end
    if any(~c1 & ~c2) % turbulent regime
        f(~c1 & ~c2) = f_work(~c1 & ~c2);
    end
end

end

function [f_work,h_f_work] = ColebrookIter(N_Re_work,f_old,epsilon,gradient)
    
help03 = 18.7./(N_Re_work.*(f_old.^(0.5)));
f_work = 1./(1.74 - 2.*(log(2.*epsilon + help03)./log(10))).^2;

if gradient
   
   % thanks adimat!
   tmp08= (0.5).* f_old.^ ((0.5)- 1);
   tmp08(f_old== 0& (0.5)== 0)= 0; % Ensure, that derivative of 0.^0 is 0 and not NaN.
   tmp00= f_old.^ (0.5);
   h_tmp01= N_Re_work.* tmp08;
   tmp01= N_Re_work.* tmp00;
   h_help03= (-18.7.* h_tmp01)./ tmp01.^ 2;
   help03= 18.7./ tmp01; 
   tmp02= 2.* epsilon;
   tmp03= tmp02+ help03;
   h_tmp_log_00000= h_help03./ tmp03;
   tmp_log_00000= log(tmp03);
   h_tmp04= h_tmp_log_00000./ log(10);
   tmp04= tmp_log_00000./ log(10);
   h_tmp05= 2.* h_tmp04;
   tmp05= 2.* tmp04;
   h_tmp06= -h_tmp05;
   tmp06= 1.74- tmp05;
   h_tmp07= 2.* tmp06.* h_tmp06;
   tmp07= tmp06.^ 2;
   h_f_work= (- h_tmp07)./ tmp07.^ 2;
    
    
end


end


