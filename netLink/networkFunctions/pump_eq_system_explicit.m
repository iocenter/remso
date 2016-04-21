function [freq ] = pump_eq_system_explicit(qf, dhf, fref, nStages, explicit)
% pump_eq_system_explicitCalculates the head difference for the pump at a given pressure loss.
% It can be used to compute upper and lower bounds for the pressure loss
% in the equipment. The coefficients of the curve were found in the paper
% 'Exploring the potential of model-based optimization in oil production
% gathering networks with esp-produced, high water cut wells'.
% A system of equations is solved to obtain the frequency of the pump for a 
% given pressure drop and flow rate flowing through it.
    if nargin < 5
        explicit = true;
    end
    
    %% attention: signH was implemented to handle negative head. If problems are found during optimization,
    %% review the equations.
    qf = -qf; % convention: positive flow for producing pump (opposite to JDJ)
    signH = sign(dhf);

    a0 = 19.37.*ones(numel(double(qf)),1);       % m 
    a1 = -1.23e-02.*ones(numel(double(qf)),1);     % in m/sm3/d
    a2 = 2.24e-05.*ones(numel(double(qf)),1);      % in m/(sm3/d)^2      
    a3 = -1.86e-08.*ones(numel(double(qf)),1);     % in m/(sm3/d)^3
    a4 = 4.13e-12.*ones(numel(double(qf)),1);      % in m/(sm3/d)^4  
    
    a2p = a2 - signH.*(dhf)./(((qf./(meter^3/day)).^2).*nStages);
    
    if explicit
        x = explicit4thPolynomial([double(a4) double(a3) double(a2p) double(a1) double(a0)]);
    else
        x = matlabRoots([double(a4) double(a3) double(a2p) double(a1) double(a0)]);
    end
    
    x = pickSolution(x, fref);
    
    [r,dxdr] = G(x,a4,a3,a2p, a1,a0);
    x = x - 1./dxdr.*r;
    freq = signH.*((qf./(meter^3/day)).*fref)./x;    
end

function y = pickSolution(x, qf)
    y = zeros(size(x,1),1);
    x = cleanImaginaryInfinitesimal(x);        
    [val, ind] = min(abs(bsxfun(@minus,x,qf)),[],2);
    
    for i=1:size(x,2);
        ik = ind ==i;
        y(ik) = x(ik,i);
    end    
end

function y = cleanImaginaryInfinitesimal(y)    
    [m,n] = size(y);
    y = reshape(y,numel(y),1);
    yImag = imag(y);
    condImag = abs(yImag)<1e-09;
    y(condImag) = real(y(condImag));
    imagNumbers = imag(y)~=0;
    y(imagNumbers) = inf;
    y = reshape(y, m,n);
end


function x = matlabRoots(coef)
    x = zeros(size(coef,1),4);
    for i=1:size(coef,1)
        x(i,:) = roots(coef(i,:));
    end
end

function [r, drdx] = G(x,c1, c2, c3, c4, c5)
%     r = c(:,1).*x.^4 + c(:,2).*x.^3 + c(:,3).*x.^2 + c(:,4).*x + c(:,5);        
    r = c1(:).*x.^4 + c2(:).*x.^3 + c3(:).*x.^2 + c4(:).*x + c5(:);    
    c = [double(c1), double(c2), double(c3), double(c4), double(c5)]; % dg/dx is cancelled in the gradient calculation    
    drdx = 4.*c(:,1).*x.^3 + 3.*c(:,2).*x.^2 + 2.*c(:,3).*x + c(:,4);
    
    assert(all(abs(real(double(r))) <= 1e-06)) 
end


function [sol] = explicit4thPolynomial(coef)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% coefficients for a 4th order polynomial eq. (quartic equation)  %%
% a4p*x^4 + a3p*x^3 + a2p*x^2 + a1p*x + aop, with x=fref/f        %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
e = coef(:,5);
d = coef(:,4);
c = coef(:,3);
b = coef(:,2);
a = coef(:,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Converting to a depressed quartic equation  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% p,q,S, and Q are terms used to solve analitically the 4th order eq %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameters to perform an analysis of the solutions  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
deltaAnalysis = 256.*(a.^3).*(e.^3) - 192.*(a.^2).*b.*d.*(e.^2) - 128.*(a.^2).*(c.^2).*(e.^2) + 144.*(a.^2).*c.*(d.^2).*e - 27.*(a.^2).*(d.^4)...
    + 144.*a.*(b.^2).*c.*(e.^2) - 6.*a.*(b.^2).*(d.^2).*e - 80.*a.*b.*(c.^2).*d.*e + 18.*a.*b.*c.*(d.^3) +  16.*a.*(c.^4).*e ...
    - 4.*a.*(c.^3).*(d.^2) - 27.*(b.^4).*(e.^2) + 18.*(b.^3).*c.*d.*e - 4.*(b.^3).*(d.^3) - 4.*(b.^2).*(c.^3).*e + (b.^2).*(c.^2).*(d.^2);

PAnalysis = 8.*a.*c - 3.*(b.^2);
QAnalysis = b.^3 + 8.*(a.^2).*d - 4.*a.*b.*c ;

% DAnalysis = 64.*(a.^3).*e - 16.*(a.^2).*(c.^2) + 16.*a.*(b.^2).*c - 16.*(a.^2).*b.*d - 3.*(b.^4);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameters to determine the solutions for the equation %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
delta0 = c.^2 - 3.*b.*d + 12.*a.*e;
delta1 = 2.*(c.^3) - 9.*b.*c.*d + 27.*(b.^2).*e + 27.*a.*(d.^2) -72.*a.*c.*e;

p = (PAnalysis)./(8.*(a.^2));
q = (QAnalysis)./(8.*(a.^3));

%     Q =  ((delta1 + (delta1.^2 -4.*(delta0.^3)).^(1/2))./(2)).^(1/3);
Q =  ((delta1 + (-27.*deltaAnalysis).^(1/2))./(2)).^(1/3);
S = 0.5.*(-(2/3).*p + (1./(3.*a)).*(Q + delta0./Q)).^(1/2);


x1 = -b./(4.*a) - S + 0.5.*(-4.*(S.^2) - 2.*p + q./S).^(1/2);
x2 = -b./(4.*a) - S - 0.5.*(-4.*(S.^2) - 2.*p + q./S).^(1/2);
x3 = -b./(4.*a) + S + 0.5.*(-4.*(S.^2) - 2.*p - q./S).^(1/2);
x4 = -b./(4.*a) + S - 0.5.*(-4.*(S.^2) - 2.*p - q./S).^(1/2);
sol = [x1, x2, x3, x4];

end
