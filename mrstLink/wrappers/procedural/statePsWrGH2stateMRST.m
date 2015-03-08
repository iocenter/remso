function [ stateMrst,Jac ] = statePsWrGH2stateMRST( p,sW,rGH,f,disgas,vapoil,varargin )
%
% Given the oil pressures (p), water saturations and gas to hidrocarbon ratio
% rGH = vG/(vG+vO) in grid blocks returns the corresponding MRST-state.
% It also returns the Jacobian of the MRST variables (pressure,sW,x)
% with respect to (p,sW,rGH) if required by the options.
% The transformation fails if 0 == (1-rvSat*rsSat) for saturated blocks.
% This condition should be tested before optimization.


opt = struct('partials',false,'tol',sqrt(eps));

opt = merge_options(opt, varargin{:});


if opt.partials
    [p, sW, rGH] = initVariablesADI(p, sW, rGH);
end


if disgas
    rsSat = f.rsSat(p);
else % otherwise rs = rsSat = const
    assert(isfield(f, 'rsSat'),'Place a constant value for the dissolved gas ratio ''rsSat'' in the fluid properties. This constant may be 0.');
    rsSat = f.rsSat;
end
if vapoil
    rvSat = f.rvSat(p);
else % otherwise rv = rvSat = const
    assert(isfield(f, 'rvSat'),'Place a constant value for the vaporized oil ratio ''rvSat'' in the fluid properties. This constant may be 0.');
    rvSat = f.rvSat;
end
rs= rsSat;
rv= rvSat;


if disgas
    [sO,sG,rsOut,st1] = undersaturatedOil(p, sW, rGH,rsSat,f.rhoOS,f.rhoGS);
    rs(st1) = rsOut(st1);
else
    % initialize variables
    sO = p*0;
    sG = p*0;
    st1 = false(size(double(p)));
end

st2 = ~st1;

if vapoil && any(st2)
    [sO(st2),sG(st2),rvOut,st2Out] = undersaturatedGas(p(st2),sW(st2),rGH(st2),rvSat(st2),f.rhoOS,f.rhoGS);
    st2(st2) = st2Out;
    rv(st2) = rvOut(st2Out);
else
    st2 = false(size(double(p)));
end

st3 = ~or(st1,st2);

if any(st3)
    if disgas
        rsIn = rsSat(st3);
    else % rsSat must be a constant
        rsIn = repmat(rsSat,sum(st3),1);
    end
    if vapoil
        rvIn = rvSat(st3);
    else % rvSat must be a constant
        rvIn = repmat(rvSat,sum(st3),1);
    end
    [sO(st3),sG(st3),st3(st3)] = saturatedFluid(p(st3),sW(st3),rGH(st3),rsIn,rvIn,f.bO,f.bG,disgas,vapoil,f.rhoOS,f.rhoGS,opt.tol);
end

assert(all(or(st1,or(st2,st3))),'Couldn''t find fluid saturation status');


stateMrst.pressure = double(p);
stateMrst.s = [double(sW),double(sO),double(sG)];
stateMrst.rs = double(rs);
stateMrst.rv = double(rv);
stateMrst.status = st1 + st2*2 + st3*3;


Jac = [];
if opt.partials
    x = st1.*rs+st2.*rv+st3.*sG;
    
    reducedState = [p;sW;x];
    Jac = reducedState.jac;
    
end


end

function [sO,sG,rs,st1] = undersaturatedOil(p, sW, rGH,rsSat,rhoOS,rhoGS)

% initialize variables
sG = p*0;   % 0 by definition
sO = 1-sW;
rs = rsSat;

st1 = ((rhoOS/rhoGS)*double(rGH) - (1-double(rGH)).*double(rsSat) < -eps);

if any(st1)
    rs(st1) = (rhoOS/rhoGS)*rGH(st1)./(1-rGH(st1));
end

end

% this function is reciprocal to undersaturatedOil
function [sO,sG,rv,st2] = undersaturatedGas(p,sW,rGH,rvSat,rhoOS,rhoGS)

% initialize variables
sG = 1-sW;
sO = p*0;  % 0 by definition
rv = rvSat;

st2 = ((rhoGS/rhoOS)*(1-double(rGH))- double(rvSat).*double(rGH) < -eps) ;

if any(st2)
    rv(st2)  =  (rhoGS/rhoOS)*(1./rGH(st2) - 1) ;
end

end




function [sO,sG,st3] = saturatedFluid(p,sW,rGH,rsSat,rvSat,bOF,bGF,disgas,vapoil,rhoOS,rhoGS,tol)
% solves  the linear system
%                       ag                                                 -ao
%[ 0    ]  =     [ bG*(rGH*(rhoGS+rhoOS*rvSat)-rhoGS), bO*(rGH(rhoOS+rhoGS*rsSat)-rhoGS*rsSat) ] * [sG]
%[ 1-sW ]        [ 1                                 , 1                                       ]   [sO]


% initialize variables
sO = p*0;


if disgas
    bO  = bOF(p, rsSat, true(size(double(p))));
else
    bO  = bOF(p);
end
if vapoil
    bG  = bGF(p, rvSat, true(size(double(p))));
else
    bG  = bGF(p);
end


ag = bG.*(rGH.*(rhoGS+rhoOS*rvSat)-rhoGS);
ao = bO.*(rhoGS*rsSat-rGH.*(rhoOS+rhoGS*rsSat));


invDetSw = (1-sW)./(ao + ag);

if any(~isfinite(double(invDetSw)))
    warning('This transformation is not appropriate for your application:  Check the transformation conditions -->  1 ~= rsSat*rvSat')
end


sG = invDetSw.*ao;
[sG,st3] = checkSaturationValueBoundary(sG,tol);  % according to theory this is a redundant check if the fluid is saturated, however numerical approximation is messing up



if any(st3)
    sO(st3) = 1-sG(st3)-sW(st3);
    [sO(st3),st3(st3)] = checkSaturationValueBoundary(sO(st3),tol);
end



end


function [s,sVF] = checkSaturationValueBoundary(s,tol)
if isa(s,'ADI')
    adi = true;
    sV = s.val;
else
    adi = false;
    sV = s;
end

sVF =  and(sV >= 0,sV <= 1);
sVT = sV(~sVF);
if any(~sVF)
   
    % ok... let's tolerate a bit if they are out of bounds
    sVFL = and(sV(~sVF) <=0,sV(~sVF) >= 0-tol);
    sVFU = and(sV(~sVF) >=1,sV(~sVF) <= 1+tol);
    
    
    if any(sVFL)
        sVT(sVFL) = 0;
    elseif any(sVFU)
        sVT(sVFU) = 1;
    end
    
    changedFlag = or(sVFL,sVFU);
    sV(~sVF) = sVT;
    sVF(~sVF)= changedFlag;
    
    if any(~sVF)
        warning('Check if fuild properties are consistent')
    end
end

if adi
    s.val = sV;
else
    s = sV;
end



end
