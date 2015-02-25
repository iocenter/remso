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


%check for p-dependent porv mult:
pvMult = 1;
if isfield(f, 'pvMultR')
    pvMult =  f.pvMultR(p);
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

V = system.s.pv.* pvMult;

if vapoil
    vO = V.* (bO.* sO  + rv.* bG.* sG);
else
    vO = V.* bO.* sO;
end

if disgas
    vG = V.*(bG.* sG  + rs.* bO.* sO);
else
    vG = V.*(bG.* sG );
end


end
