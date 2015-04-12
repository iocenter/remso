function  [alpha,f,g,vars,simVars,neval,xfd,debugInfo] = linesearch(fun,f0,g0,f1,g1,eta,varargin)
% Line-search
%
% SYNOPSIS:
%  [alpha,f,g,vars,simVars,neval,xfd,debugInfo] = linesearch(fun,f0,g0,f1,g1,eta)
%  [alpha,f,g,vars,simVars,neval,xfd,debugInfo] = linesearch(fun,f0,g0,f1,g1,eta, 'pn', pv, ...)
%
% PARAMETERS:
%
%   fun - line function
%
%   (f0,g0) - function and gradient evaluted on 0
%
%   (f1,g1) - function and gradient evaluated on 1
%
%   eta -  sufficient decrease parameter
%
%   'pn'/pv - List of 'key'/value pairs defining optional parameters. The
%             supported options are:
%
%   tau - parameter ralated to the gradient condition on the line-search
%
%   kmax - maximum number of extra points to be evaluated
%
%   curvLS - if true then honor the curvature condition during line-search
%
%   required - required decrease of the merit function
%
%   skipWatchDog - Perform normal line-search.
%
%   debug - true to collect debug information
%
%
% RETURNS:
%
%
%   (alpha,f,g,) - Line value, function value and gradient obtained
%
%  (vars,simVars) - inputs and outputs of fun at the line-search solution
%  that are returned to avoid recomputation
%
%   neval - number of point evaluations during line-search
%
%   xfd - record of line value, function, and gradients during line-search
%
%   debugInfo - debug information
%
% SEE ALSO:
%
%


opt = struct('tau',0.1,'kmax',3,'debug',false,'curvLS',true,'required',inf,'maxStep',1);
opt = merge_options(opt, varargin{:});

fvalb = []; %% expected value by polynomial model

% Initial values
alpha = -1;
f=0;
g=0;
neval = 0;
xfd = zeros(opt.kmax,3);
vars = [];
simVars = [];
debugInfo = cell(opt.kmax,1);


% Check descent condition
if  g0 >= 0
    xfd = zeros(0,3);
    debugInfo = {};
    warning('linesearch initialized without a decress direction')
    return
end

% Conditions set-up
armijoF = @(lT,fT)  (fT - (f0 + eta*g0*lT));
armijoOk = @(lT,fT) (armijoF(lT,fT) <= 0);
if opt.curvLS
    curvatureOk = @(gT) (gT >= opt.tau*g0);
else
    curvatureOk = @(gT) true;
end


%Armijo is not ok! otherwise we wouldn't be here!
stepb = opt.maxStep;        
fib = f1;
gb = g1;

xfdS = [0,f0,g0;
        stepb,f1,g1];
while  (neval < opt.kmax) && ( ~curvatureOk(gb) ||  ~armijoOk(stepb,fib) || (f > opt.required) ) 
    
    
    [stepb,fvalb,inMin,xfdS] = interpolate(xfdS);
    if ~inMin || fvalb > opt.required  %% if you believe that you can make it, continue
        break;
    end
    
    [fib , gb,v,sV,debugInfoN] = fun(stepb);
    
    neval = neval+1;
    xfd(neval,:) = [stepb  fib  gb];
    debugInfo{neval} = debugInfoN;
    debugInfo{neval}.armijoVal = armijoF(stepb,fib);
    
    if opt.debug
        plotData(f0,g0,f1,g1,xfd,neval)
        plot(stepb,fvalb,'xr');
    end
    
    xfdS = [xfdS;stepb,fib,gb];
    if ( armijoOk(stepb,fib) )
        % move incumbet step
        alpha = stepb;
        f = fib;
        g = gb;
        vars = v;
        simVars = sV;
    end
    
end

if opt.debug
    plotData(f0,g0,f1,g1,xfd,neval)
    if ~isempty(fvalb)
        plot(stepb,fvalb,'xr');
    end
end

% Return values

xfd = xfd(1:neval,:);
debugInfo = debugInfo(1:neval);

end



function  [alpha,fval,inMin,xfdS] = interpolate(xfd)
%  find minimizer of data



inMin = false;
n = size(xfd,1);
if n <= 1  % so... you want to interpolate just giving one point?
    error('There are no enough points to construct model')
end
if n == 2
    xfdS = sortrows(xfd(1:2,:));
    coef = constructCubicPolynomial(xfdS(1:2,:));
    [alpha,fval,inMin] = findMinimumCubicPol(xfdS(1,1),xfdS(2,1),coef);
else  %% n == 3
    xfdS = sortrows(xfd([n-1,n],:));
    coef = constructCubicPolynomial(xfdS(1:2,:));
    [alpha1,fval1,inMin1] = findMinimumCubicPol(xfdS(1,1),xfdS(2,1),coef);
    
    xfdS = sortrows(xfd([n-2,n],:));
    coef = constructCubicPolynomial(xfdS(1:2,:));
    [alpha2,fval2,inMin2] = findMinimumCubicPol(xfdS(1,1),xfdS(2,1),coef);
    
    if fval1 < fval2
        alpha = alpha1;
        fval = fval1;
        inMin = inMin1;
    else
        alpha = alpha2;
        fval = fval2;
        inMin = inMin2;
    end
    
    
end





end


function coef = constructCubicPolynomial(xfd2)
%
% construct a cubic Polynomial fitting to xfd2
%

C = zeros(4,4);
d = zeros(4,1);


for i = 1:2
    
    C(2*i-1,1) =   1;
    C(2*i,1)   =   0;
    for j = 2:4
        C(2*i-1,j) = xfd2(i,1)^(j-1);
        C(2*i,j) =  (j-1)*xfd2(i,1)^(j-2);
    end
    
    d(2*i-1)  = xfd2(i,2);
    d(2*i)  = xfd2(i,3);
end


coef = lsqlin(C,d);

%{

pol = @(x) sum(coef'.*x.^[0,1,2,3]);

x = (0:1:100) * (xfd2(2,1)-xfd2(1,1))/100+xfd2(1,1);

y = zeros(size(x));
for k = 1:numel(x)
    y(k) = pol(x(k))
end
plot(x,y)

%}


end

function [a,fval,inMin] = findMinimumCubicPol(xm,xM,coef)
%
% find the minimum of a third-order polynomial on [xm xM]
%

inMin = false;

pol = @(x) sum(coef'.*x.^[0,1,2,3]);

vals = [pol(xm),xm;
    pol(xM),xM];

[fval,ind] = min(vals(:,1));
a = vals(ind,2);

if coef(4) ~= 0 %% polynomial is cubic find both roots
    
    gradCoef = [coef(2),2*coef(3),3*coef(4)];
    
    rs = gradCoef(2)^2-4*gradCoef(3)*gradCoef(1);
    
    if rs >= 0  % no real roots
        signs = [1,-1];
        
        for s = signs
            xr = (-gradCoef(2)+s*sqrt(rs))/(2*gradCoef(3));
            if and((xm < xr),(xr < xM))
                polval = pol(xr);
                if polval < fval
                    inMin =true;
                    fval = polval;
                    a = xr;
                end
            end
        end
    end
    
else %% polynomial is quadratic
    if coef(3) ~= 0  %% polynomial is quadratic
        xr = coef(2)/(2*coef(3));
        if and((xm < xr),(xr < xM))
            polval = pol(xr);
            if polval < fval
                inMin =true;
                fval = polval;
                a = xr;
            end
        end
    else  % if it is linear there is nothing to do!
        
    end
    
end


end


function [] = plotData(f0,g0,f1,g1,xfd,n)

xfd = [0,f0,g0;
       1,f1,g1;
       xfd];

n= n +2;


figure(20)
clf;
hold on;
xmax = max(xfd(:,1));
slopex = @(x) [x-xmax/10,x+xmax/10];
slopef = @(f,g) [f-g*xmax/10,f+g*xmax/10];

plot(xfd(1:n,1),xfd(1:n,2),'x')

for k = 1:n
    plot(slopex(xfd(k,1)),slopef(xfd(k,2),xfd(k,3)))
end



end


