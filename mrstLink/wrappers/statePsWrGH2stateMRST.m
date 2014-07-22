function [ stateMrst,Jac ] = statePsWrGH2stateMRST( p,sW,rGH,f,system,varargin )
%
% Given the oil pressures (p), water saturations and gas to hidrocarbon ratio
% rGH = vG/(vG+vO) in grid blocks returns the corresponding MRST-state.
% It also returns the Jacobian of the MRST variables (pressure,sW,x)
% with respect to (p,sW,rGH) if required by the options.
% The transformation fails if ....


opt = struct('partials',false);

opt = merge_options(opt, varargin{:});

comp = system.activeComponents;

disgas= comp.disgas;
vapoil= comp.vapoil;

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
    [sO,sG,rsOut,st1] = undersaturatedOil(p, sW, rGH,rsSat);
    rs(st1) = rsOut(st1);
else
    % initialize variables
    sO = p*0;
    sG = p*0;
    st1 = false(size(double(p)));
end

st2 = ~st1;

if vapoil && any(st2)
    [sO(st2),sG(st2),rvOut,st2Out] = undersaturatedGas(p(st2),sW(st2),rGH(st2),rvSat(st2));
    rv(st2) = rvOut(st2Out);
    st2(st2) = st2Out;
else
    st2 = false(size(double(p)));
end

st3 = ~or(st1,st2);

if any(st3)
    if disgas
        rsIn = rsSat(st3);
    else % rsSat must be a constant
        rsIn = ones(sum(st3),1).*rsSat;
    end
    if vapoil
        rvIn = rvSat(st3);
    else % rvSat must be a constant
        rvIn = ones(sum(st3),1).*rvSat;
    end
    [sO(st3),sG(st3),st3(st3)] = saturatedFluid(p(st3),sW(st3),rGH(st3),rsIn,rvIn,f.bO,f.bG,disgas,vapoil);
end

if ~all(or(st1,or(st2,st3)))
    assert(all(or(st1,or(st2,st3))),'Couldn''t find fluid saturation status');
end

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

function [sO,sG,rs,st1] = undersaturatedOil(p, sW, rGH,rsSat)

% initialize variables
sG = p*0;   % 0 by definition
sO = 1-sW;
rs = rsSat;

st1 = (rGH <= (1-rGH).*rsSat);

if any(st1)
    rs(st1) = rGH(st1)./(1-rGH(st1));
end

end

% this function is reciprocal to undersaturatedOil
function [sO,sG,rv,st2] = undersaturatedGas(p,sW,rGH,rvSat)

% initialize variables
sG = 1-sW;
sO = p*0;  % 0 by definition
rv = rvSat;

st2 = (1-rGH <= rvSat.*rGH) ;

if any(st2)
    rv(st2)  =  1./rGH(st2) - 1 ;
end

end




function [sO,sG,st3] = saturatedFluid(p,sW,rGH,rsSat,rvSat,bOF,bGF,disgas,vapoil)
% solves  the linear system
%                       ag                         -ao
%[ 0    ]  =     [ bG*(rGH*(1+rvSat)-1), bO*(rGH(1+rsSat)-rsSat) ] * [sG]
%[ 1-sW ]        [ 1                   , 1                       ]   [sO]


% initialize variables
sO = p*0;

% TODO: Call saturated function after bug fix.
if disgas
    bO  = bOF(p, rsSat, false(size(double(p))));
else
    bO  = bOF(p);
end
if vapoil
    bG  = bGF(p, rvSat, false(size(double(p))));
else
    bG  = bGF(p);
end


ag = bG.*(rGH.*(1+rvSat)-1);
ao = bO.*(rsSat-rGH.*(1+rsSat));


invDetSw = (1-sW)./(ao + ag);

if any(~isfinite(double(invDetSw)))
    warning('This transformation is not appropriate for your application:  Check the transformation conditions -->  1 ~= rsSat*rvSat')
end


sG = invDetSw.*ao;
st3 = and(sG >= 0,sG <= 1);  % according to theory this is a redundant check if the fluid is saturated

if any(st3)
    sO(st3) = 1-sG(st3)-sW(st3);
    st3(st3) = and(sO(st3) >= 0,sO(st3) <= 1); % according to theory this is a redundant check if the fluid is saturated
end



end
