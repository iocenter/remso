function [ stateMrst,Jac ] = statePOG2stateMRST( p,vO,vG,f,system,varargin )
%
% Given the oil pressures (p) in grid blocks and the std volumes of oil and gas
% (vO,vG) returns the corresponding MRST-state.  It also returns the
% Jacobian of the MRST variables (pressure,sW,x) with respect to (p,vO,vG)
% if required by the options.
% The transformation fails if 0 == (1-rvSat*rsSat) for saturated blocks.
% This condition should be tested before optimization.
opt = struct('partials',false);

opt = merge_options(opt, varargin{:});

comp = system.activeComponents;

disgas= comp.disgas;
vapoil= comp.vapoil;

if opt.partials
    [p, vO, vG] = initVariablesADI(p, vO, vG);
end

%check for p-dependent porv mult:
pvMult = 1;
if isfield(f, 'pvMultR')
    pvMult =  f.pvMultR(p);
end
V = system.s.pv.* pvMult;


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
    [sW,sO,sG,rsOut,st1] = undersaturatedOil(p,vO,vG,V,rsSat,f.bO);
    rs(st1) = rsOut(st1);
else
    % initialize variables
    sW = p*0;
    sO = p*0;
    sG = p*0;
    rs = rsSat;
    st1 = false(size(double(p)));
end

st2 = ~st1;

if vapoil && any(st2)
    [sW(st2),sO(st2),sG(st2),rvOut,st2Out] = undersaturatedGas(p(st2),vO(st2),vG(st2),V(st2),rvSat(st2),f.bG);
    rv(st2) = rvOut(st2Out);
    st2(st2) = st2Out;
else
    st2 = false(size(double(p)));
end

st3 = ~or(st1,st2);

if any(st3)
    if disgas
        rsIn = rsSat(st3);
    else
        rsIn = ones(sum(st3),1).*rsSat;
    end
    if vapoil
        rvIn = rvSat(st3);
    else
        rvIn = ones(sum(st3),1).*rvSat;
    end
    [sW(st3),sO(st3),sG(st3),st3(st3)] = saturatedFluid(p(st3),vO(st3),vG(st3),V(st3),rsIn,rvIn,f,disgas,vapoil);   
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

function [sW,sO,sG,rs,st1] = undersaturatedOil(p,vO,vG,V,rsSat,b)

% initialize variables
sG = p*0;   % 0 by definition
sO = p*0;
sW = p*0;
rs = rsSat;
st1 = false(size(double(p)));

oilFlag = (vO>0);

rs(oilFlag) = vG(oilFlag)./vO(oilFlag);

st1(oilFlag) = (rs(oilFlag)<rsSat(oilFlag));

bO  = b(p(st1), rs(st1), ~st1(st1));  %% ~st1(st1) == false
sO(st1) = vO(st1)./(V(st1).*bO);

st1(st1) = (sO(st1) <= 1);

sW(st1) = 1 - sO(st1);

% This commented function is the same as above but not vectorized!
%sG = 0;
%rv = rvSat;
%
% if vO > 0
%     rs = vG/vO;
%     st1 = (rs < rsSat);
%     if st1
%         bO  = b(p, rs, false(size(double(p))));
%         sO = vO/(V*bO);
%         sW = 1 - sO;
%     else
%         % THIS IS NOT UNDERSATURATED!
%         sO = 0;
%         sW = 0;
%     end
%
% else  %% vO == 0
%     st1 = false;
%     sO = 0;
%     sW = 1;
%     rs = rsSat;
% end

end

% this function is reciprocal to undersaturatedOil
function [sW,sO,sG,rv,st2] = undersaturatedGas(p,vO,vG,V,rvSat,b)

% initialize variables
sG = p*0;
sO = p*0;  % 0 by definition
sW = p*0;
rv = rvSat;
st2 = false(size(double(p)));

gasFlag = (vG>0);

rv(gasFlag) = vO(gasFlag)./vG(gasFlag);

st2(gasFlag) = (rv(gasFlag)<rvSat(gasFlag));

bG  = b(p(st2), rv(st2), ~st2(st2));   %% ~st2(st2) == false
sG(st2) = vG(st2)./(V(st2).*bG);

st2(st2) = (sG(st2) <= 1);

sW(st2) = 1 - sG(st2);


end




function [sW,sO,sG,st3] = saturatedFluid(p,vO,vG,V,rsSat,rvSat,f,disgas,vapoil)
% solves  the linear system
%[ vO ]  = V*  [ bo       , rvSat*bg ] * [sO]
%[ vG ]        [ rsSat*bo ,  bg      ]   [sG]


% initialize variables
sO = p*0;
sW = p*0;

% TODO: Call saturated function after bug fix.
if disgas
    bO  = f.bO(p, rsSat, false(size(double(p))));
else
    bO  = f.bO(p);
end
if vapoil
    bG  = f.bG(p, rvSat, false(size(double(p))));
else
    bG  = f.bG(p);
end


invDet = (1./(V.*bO.*bG.*(1-rsSat.*rvSat)));

if any(~isfinite(double(invDet)))
    warning('This transformation is not appropriate for your application:  Check that rsSat*rvSat ~= 1 for saturated fluids!')
end


sG = invDet.*(-rsSat.*bO    .*vO  +bO    .*vG);
st3 = and(sG >= 0,sG <= 1);

sO(st3) = invDet(st3).*(bG(st3)  .*vO(st3)  -rvSat(st3).*bG(st3)   .*vG(st3));
st3(st3) = and(sO(st3) >= 0,sO(st3) <= 1);


sW(st3) = 1 - sO(st3) - sG(st3);
st3(st3) = and(sW(st3) >= 0,sW(st3) <= 1);



% if st3
%     sO = invDet*(bg  *vO  -rv*bG   *vG);
%     sW = 1 - sO - sG;
% else
%     % fluid is not saturated!
%     sO = 0;
%     sW = 0;
% end




end
