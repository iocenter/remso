function [ p,rO,rG ] = stateMrst2statePOG(state,f,disgas,vapoil,varargin)
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
[st1,st2,st3]  = getCellStatus(state,  disgas, vapoil);

% define primary varible x and initialize
x = st1.*rs + st2.*rv + st3.*sG;

if opt.partials
    [p, sW, x] = initVariablesADI(p, sW, x);
end

% calculateHydrocarbonsFromStatus
% define sG, rs and rv in terms of x
sG = st2.*(1-sW) + st3.*x;
if disgas
    rsSat = f.rsSat(p);
    rs = (~st1).*rsSat + st1.*x;
else % otherwise rs = rsSat = const
    %rsSat = rs;
end
if vapoil
    rvSat = f.rvSat(p);
    rv = (~st2).*rvSat + st2.*x;
else % otherwise rv = rvSat = const
    %rvSat = rv;
end

% OIL PROPS
if disgas
    %bO  = f.bO(p, rs, false(size(double(p))));?
	bO  = f.bO(p,  rs, ~st1);
else
    bO  = f.bO(p);
end


% GAS PROPS (calculated at oil pressure)
if vapoil
    %bG  = f.bG(p, rv, false(size(double(p))));?
	bG  = f.bG(p, rv, ~st2);
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


end

function [st1, st2, st3] = getCellStatus(state, disgas, vapoil)
% Status should be passed on from updateStateVO (to be sure definition is
% identical). rs and rv are assumed to be compatible, i.e. rx = rxSat for
% saturated cells and rx <= rxSat for undersaturated. Three values of
% status are:
% status 0: should not occur (almost water only -> state 3)
% status 1 oil, no gas  : x = rs, sg = 0    , rv = rvMax
% status 2 gas, no oil  : x = rv, sg = 1-sw , rs = rsMax
% status 3 oil and gas  : x = sg, rs = rsMax, rv = rvMax
if isfield(state, 'status')
    status = state.status;
else
    s = state.s;
    watOnly    = s(:,1) > 1- sqrt(eps);
    if ~vapoil
        oilPresent = true;
    else
        oilPresent = or(s(:,2) > 0, watOnly);
    end
    if ~disgas
        gasPresent = true;
    else
        gasPresent = or(s(:,3) > 0, watOnly);
    end
    status = oilPresent + 2*gasPresent;
end
if ~disgas
    st1 = false;
else
    st1 = status==1;
end
if ~vapoil
    st2 = false;
else
    st2 = status==2;
end
st3 = status == 3;
end
