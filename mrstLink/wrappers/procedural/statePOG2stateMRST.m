function [ stateMrst,Jac ] = statePOG2stateMRST( p,rO,rG,f,disgas,vapoil,varargin )
%
% Given the oil pressures (p) in grid blocks and the adimensional volumes of oil and gas
% (rO,rG) returns the corresponding MRST-state.  It also returns the
% Jacobian of the MRST variables (pressure,sW,x) with respect to (p,rO,rG)
% if required by the options.
% The transformation fails if 0 == (1-rvSat*rsSat) for saturated blocks.
% This condition should be tested before optimization.


opt = struct('partials',false);

opt = merge_options(opt, varargin{:});



if opt.partials
    [p, rO, rG] = initVariablesADI(p, rO, rG);
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
    [sW,sO,sG,rsOut,st1] = undersaturatedOil(p,rO,rG,rsSat,f.bO,f.rhoOS,f.rhoGS);
    rs(st1) = rsOut(st1);
else
    % initialize variables
    sW = p*0;
    sO = p*0;
    sG = p*0;
    st1 = false(size(double(p)));
end

st2 = ~st1;

if vapoil && any(st2)
    [sW(st2),sO(st2),sG(st2),rvOut,st2Out] = undersaturatedGas(p(st2),rO(st2),rG(st2),rvSat(st2),f.rhoOS,f.rhoGS);
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
    [sW(st3),sO(st3),sG(st3),st3(st3)] = saturatedFluid(p(st3),rO(st3),rG(st3),rsIn,rvIn,f.bO,f.bG,disgas,vapoil,f.rhoOS,f.rhoGS);   
end

%{
st4 = ~or(st1,or(st2,st3));  %unclassified

% This transformation works, but it is difficult to bound the set of values 
p,ro,rg that have a corrsponding solution in the mrst states!

explain the paper!

[sWU,sOU,sGU,rsOutU,st1U] = undersaturatedOil(p(st4),rO(st4),rG(st4),rsSat(st4),f.bO,f.rhoOS,f.rhoGS);



    if disgas
        rsIn = rsSat(st4);
    else % rsSat must be a constant
        rsIn = repmat(rsSat,sum(st3),1);
    end
    if vapoil
        rvIn = rvSat(st4);
    else % rvSat must be a constant
        rvIn = repmat(rvSat,sum(st4),1);
    end

[sWU,sOU,sGU,st3U = saturatedFluid(p(st4),rO(st4),rG(st4),rsIn,rvIn,f.bO,f.bG,disgas,vapoil,f.rhoOS,f.rhoGS);   


%}

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

function [sW,sO,sG,rs,st1] = undersaturatedOil(p,rO,rG,rsSat,b,rhoOS,rhoGS)

% initialize variables
sG = p*0;   % 0 by definition
sO = p*0;
sW = p*0;
rs = rsSat;
st1 = false(size(double(p)));

oilFlag = (rO>0);

rs(oilFlag) = (rhoOS/rhoGS)*rG(oilFlag)./rO(oilFlag);

st1(oilFlag) = (rs(oilFlag)<rsSat(oilFlag));

bO  = b(p(st1), rs(st1), ~st1(st1));  %% ~st1(st1) == false
sO(st1) = rO(st1)./bO/rhoOS;

st1(st1) = (sO(st1) <= 1);

if any(st1)
    sW(st1) = 1 - sO(st1);
end

end

% this function is reciprocal to undersaturatedOil
function [sW,sO,sG,rv,st2] = undersaturatedGas(p,rO,rG,rvSat,b,rhoOS,rhoGS)

% initialize variables
sG = p*0;
sO = p*0;  % 0 by definition
sW = p*0;
rv = rvSat;
st2 = false(size(double(p)));

gasFlag = (rG>0);

rv(gasFlag) = (rhoGS/rhoOS)*rO(gasFlag)./rG(gasFlag);

st2(gasFlag) = (rv(gasFlag)<rvSat(gasFlag));

bG  = b(p(st2), rv(st2), ~st2(st2));   %% ~st2(st2) == false
sG(st2) = rG(st2)./(bG)/rhoGS;

st2(st2) = (sG(st2) <= 1);

if any(st2)
    sW(st2) = 1 - sG(st2);
end

end




function [sW,sO,sG,st3] = saturatedFluid(p,rO,rG,rsSat,rvSat,bOF,bGF,disgas,vapoil,rhoOS,rhoGS)
% solves  the linear system
%[ rO ]  =     [ bo*rhoOS       , rvSat*bg*rhoOS ] * [sO]
%[ rG ]        [ rsSat*bo*rhoGS ,  bg*rhoGS      ]   [sG]


% initialize variables
sO = p*0;
sW = p*0;

% TODO: Call saturated function after bug fix.
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


invDet = (1./(rhoOS*rhoGS*bO.*bG.*(1-rsSat.*rvSat)));

if any(~isfinite(double(invDet)))
    warning('This transformation is not appropriate for your application:  Check that rsSat*rvSat ~= 1 for saturated fluids!')
end


sG = invDet.*(-rsSat.*bO*rhoGS    .*rO  +bO*rhoOS    .*rG);
st3 = and(sG >= 0,sG <= 1);

sO(st3) = invDet(st3).*(bG(st3)*rhoGS  .*rO(st3)  -rvSat(st3).*bG(st3)*rhoOS   .*rG(st3));
st3(st3) = and(sO(st3) >= 0,sO(st3) <= 1);


sW(st3) = 1 - sO(st3) - sG(st3);
st3(st3) = and(sW(st3) >= 0,sW(st3) <= 1);








end
