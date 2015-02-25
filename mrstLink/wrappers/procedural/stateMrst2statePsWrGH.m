function [ p,sW,rGH ] = stateMrst2statePsWrGH(state,f,disgas,vapoil,varargin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

opt = struct('partials',false);

opt = merge_options(opt, varargin{:});


% current variables: ------------------------------------------------------
p    = state.pressure;
sW   = state.s(:,1);
sG   = state.s(:,3);
rs   = state.rs;
rv   = state.rv;

%Initialization of primary variables ----------------------------------
st  = getCellStatusVO(state,  1-sW-sG,   sW,  sG,  disgas, vapoil);

% define primary varible x and initialize
x = st{1}.*rs + st{2}.*rv + st{3}.*sG;

if opt.partials
    [p, sW, x] = initVariablesADI(p, sW, x);
end

% calculateHydrocarbonsFromStatus
% define sG, rs and rv in terms of x
sG = st{2}.*(1-sW) + st{3}.*x;
if disgas
    rsSat = f.rsSat(p);
    rs = (~st{1}).*rsSat + st{1}.*x;
else % otherwise rs = rsSat = const
    %rsSat = rs;
end
if vapoil
    rvSat = f.rvSat(p);
    rv = (~st{2}).*rvSat + st{2}.*x;
else % otherwise rv = rvSat = const
    %rvSat = rv;
end


% OIL PROPS
if disgas
    %bO  = f.bO(p, rs, false(size(double(p))));?
	bO  = f.bO(p,  rs, ~st{1});
else
    bO  = f.bO(p);
end

% GAS PROPS (calculated at oil pressure)
if vapoil
    %bG  = f.bG(p, rv, false(size(double(p))));?
	bG  = f.bG(p, rv, ~st{2});

else
    bG  = f.bG(p);
end

% EQUATIONS -----------------------------------------------------------
sO  = 1- sW  - sG;


if vapoil
    rO = (bO.* sO  + rv.* bG.* sG)*f.rhoOS;
else
    rO = (bO.* sO)*f.rhoOS;
end

if disgas
    rG = (bG.* sG  + rs.* bO.* sO)*f.rhoGS;
else
    rG = (bG.* sG)*f.rhoGS;
end

rGH = p*0; % initialization

rH = rO+rG;

st4 = rH > 0;

rGH(st4) = rG(st4)./rH(st4);


% This is a remedy to allow rs=rsSat or rv=rvSat when rO+rG==0 (and therefore sW = 1):
% we take the rGH = lim_{rO -> 0} rsSat*rO/(rsSat*rO+rO) = rsSat/(1+rsSat)
if any(~st4)
    if disgas
        rGH(~st4) = f.rhoGS*rsSat(~st4)./(f.rhoOS+f.rhoGS*rsSat(~st4));
    elseif vapoil
        rGH(~st4) = f.rhoGS./(f.rhoGS+f.rhoOS*rvSat(~st4));
    else
        if f.rsSat > 0
            rGH(~st4) = rGH(~st4) +  f.rhoGS*f.rsSat./(f.rhoOS+f.rhoGS*f.rsSat);
        else
            rGH(~st4) = rGH(~st4) + f.rhoOS./(f.rhoGS+f.rhoOS*f.rvSat);
        end
    end
end




end
