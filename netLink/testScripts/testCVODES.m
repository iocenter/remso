clear; clc;

% Required MRST modules
mrstModule add deckformat
mrstModule add ad-fi

addpath(genpath('../../optimization/remso'));
addpath(genpath('../../optimization/utils'));
addpath(genpath('../../mrstDerivated'));
addpath(genpath('../../netLink'));
addpath(genpath('../../netLink/graphFunctions'));
addpath(genpath('../../netLink/networkFunctions'));
addpath(genpath('../../netLink/fluidProperties'));
addpath(genpath('../../netLink/pipeFlow'));
addpath(genpath('../../netLink/conversionFactors'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% Edge e1 %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

v1 = newVertex(1, -1,-1);
v2  = newVertex(2, -1, -1);
% v2.pressure = 800*psia; % in Pa

% v2.pressure = 720*psia; % in Pa

v2.pressure = 350*barsa;

e1 = newEdge(1, v1, v2, -1);
e1.units = 0; % METRIC=0, FIELD = 1,
% e1.pipeline = newPipeline('diam', 2.5*inch, 'len', 500*meter , 'ang', degtorad(90), 'temp',  convtemp(175,'F','K'));

e1.pipeline = newPipeline('diam', 0.249*ft, 'len', 1*ft , 'ang', degtorad(90), 'temp',  convtemp(175,'F','K'));

e1.stream = newStream();
% e1.stream.gas_visc = 0.0131*(centi*poise); % viscosity in kilogram/(meter*second)
e1.stream.gas_visc = 0.018*(centi*poise); % viscosity in kilogram/(meter*second)
e1.stream.sg_gas = 0.65;
% e1.stream.gas_dens = 2.6*pound/ft^3; % in kg/m^3
e1.stream.gas_dens = 0.8; % in kg/m^3

% e1.stream.oil_visc = 2*(centi*poise); % viscosity in kilogram/(meter*second)
e1.stream.oil_visc = 18*(centi*poise); % viscosity in kilogram/(meter*second)
e1.stream.sg_oil = 0.35;
% e1.stream.oil_dens =  49.9*pound/ft^3; % in kg/m^3
e1.stream.oil_dens =  56.6*pound/ft^3; % in kg/m^3

f = 519.68;
wcut = 0.15;

e1.qoE = -(1-wcut)*f*(meter^3/day);   % sm3/s
e1.qwE = -wcut*f*(meter^3/day);       % sm3/s
e1.qgE = -1583.84*(meter^3/day)*0;             % sm3/s


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% Edge e2 %%%%%%%%%%%%%%%%%%%%%f
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

v3 = newVertex(3, -1,-1);
v4  = newVertex(4, -1, -1);
v4.pressure = 800*psia; % in Pa

e2 = newEdge(1, v3, v4, -1);
e2.units = 0; % METRIC=0, FIELD = 1,
e2.pipeline = e1.pipeline;
e2.pipeline.diam = 2.5*inch;

e2.stream = e1.stream;
e2.stream.oil_dens = 49.9*(pound/ft^3);
e2.stream.oil_visc = 2*(centi*poise);
e2.stream.gas_dens = 2.6*(pound/ft^3);
e2.stream.oil_visc = 0.0131*(centi*poise);


% f = 2000;
% wcut = 0.0;

% e2.qoE = -(1-wcut)*f*(meter^3/day);   % sm3/s
e2.qoE = -2000*(stb/day);
e2.qwE = -250*(stb/day);
% e2.qwE = -wcut*f*(meter^3/day);       % sm3/s
% e2.qgE = -10^3*(ft^3/day);             % sm3/s
e2.qgE = -1583.84*(meter^3/day);

% E = [e1; e2];
% V = [v1; v2; v3; v4];

E = [e1];
V = [v1; v2];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Graph  0described as vectors of flows (qo, qw, qg) and pressures (p)  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Vout  = [];
for i=1:length(V)
    if V(i).id == 2 || V(i).id ==4
        Vout = [Vout; V(i)];
    end
end

[qo, qw, qg, p] = graph2FlowsAndPressures(Vout, E);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Calculating pressure drops %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dp = dpBeggsBrillJDJ(E, qo, qw, qg, p);
dp/barsa;

E.integrationStep = 0.001;

[errorMax, errorJo, errorJw, errorJg, errorJp] = testRunNetworkADI(E, qo, qw, qg, p, 'qopert', 1e-03*meter^3/day, 'qwpert', 1e-03*meter^3/day, 'qgpert',  0*meter^3/day, 'poutPert', 1e-03*barsa);


tic
[ dpTotalCVODES ] = dpCVODES(E,  qo, qw, qg, p,'dpFunction', @dpBeggsBrillJDJ)
toc

[qo, qw, qg,p] = initVariablesADI(qo, qw, qg,p);

tic
[ dpTotalCVODESFwd ] = dpCVODES(E,  qo, qw, qg, p,'dpFunction', @dpBeggsBrillJDJ)
toc

tic
[ dpTotalCVODESBack ] = dpCVODES(E,  qo, qw, qg, p,'forwardGradient',false,'dpFunction', @dpBeggsBrillJDJ)
toc


tic
[ dpTotalCVODESFwdFD ] = dpCVODES(E,  qo, qw, qg, p,'forwardGradient',true, 'finiteDiff', true, 'dpFunction', @dpBeggsBrillJDJ, 'hasSurfaceGas', false)
toc

tic
[ dpTotalCVODESBackFD ] = dpCVODES(E,  qo, qw, qg, p,'forwardGradient', false, 'finiteDiff', true, 'dpFunction', @dpBeggsBrillJDJ, 'hasSurfaceGas', false)
toc

tic
[ dpTotalStepwise ] = dpStepwise(E,  qo, qw, qg, p,'dpFunction', @dpBeggsBrillJDJ)
toc

double(dpTotalCVODES-dpTotalCVODESFwd)
dpTotalCVODESFwd-dpTotalCVODESBack
dpTotalCVODESFwd-dpTotalStepwise


dpF = @(Ei,  qoi, qwi, qgi, pi) dpCVODES(Ei,  qoi, qwi, qgi, pi,'forwardGradient',false,'dpFunction', @dpBeggsBrillJDJ);

[errorMax, joRelError, jwRelError, jgRelError, jpReError] = testRunNetworkADI(E, double(qo), double(qw), double(qg), double(p), ...
    'qopert', 1e-03*meter^3/day, 'qwpert', 1e-03*meter^3/day, 'qgpert',  0*1e-03*meter^3/day, 'poutPert', 1e-03*barsa,...
    'dpFunction',dpF)

qoBounds = [0;1000*(meter^3/day)];
qwBounds = [0;1000*(meter^3/day)];
qgBounds = [0;10000*(meter^3/day)];  
pBounds = [5*barsa;250*barsa];  

mask = [true;true;true;true];

nGrid = 10;
qoRandom = qoBounds(1) + (qoBounds(2)-qoBounds(1))*rand(nGrid,1);
qwRandom = qwBounds(1) + (qwBounds(2)-qwBounds(1))*rand(nGrid,1);
qgRandom = qgBounds(1) + (qgBounds(2)-qgBounds(1))*rand(nGrid,1);
pRandom =  pBounds(1) + ( pBounds(2)- pBounds(1))*rand(nGrid,1);




pert = 1e-5;
qoP = 5*meter^3/day*pert;
qwP = 5*meter^3/day*pert;
qgP = (10*ft)^3/day*pert;
pP =  5*barsa*pert;

gradScale = 5*barsa./[5*meter^3/day,5*meter^3/day,100*(10*ft)^3/day,5*barsa];

dpErrAbs1 = zeros(nGrid,1);
for k = 1:nGrid
    [qok, qwk, qgk,pk] = deal(qoRandom(k), qwRandom(k), qgRandom(k), pRandom(k));
    
    dpGradNum =[
        (dpCVODES(E, qok+qoP, qwk, qgk,pk , 'dpFunction', @dpBeggsBrillJDJ)-dpCVODES(E, qok-qoP, qwk, qgk,pk , 'dpFunction', @dpBeggsBrillJDJ))/(2*qoP),...
        (dpCVODES(E, qok, qwk+qwP, qgk,pk , 'dpFunction', @dpBeggsBrillJDJ)-dpCVODES(E, qok, qwk-qwP, qgk,pk , 'dpFunction', @dpBeggsBrillJDJ))/(2*qwP),...
        (dpCVODES(E, qok, qwk, qgk+qgP,pk , 'dpFunction', @dpBeggsBrillJDJ)-dpCVODES(E, qok, qwk, qgk-qgP,pk , 'dpFunction', @dpBeggsBrillJDJ))/(2*qgP),...
        (dpCVODES(E, qok, qwk, qgk,pk+ pP , 'dpFunction', @dpBeggsBrillJDJ)-dpCVODES(E, qok, qwk, qgk,pk- pP , 'dpFunction', @dpBeggsBrillJDJ))/(2* pP),...
        ];
    
    
    [qok, qwk, qgk,pk] = initVariablesADI(qoRandom(k), qwRandom(k), qgRandom(k), pRandom(k));
    dpk = dpCVODES(E, qok, qwk, qgk,pk , 'dpFunction', @dpBeggsBrillJDJ);
    dpGrad = cell2mat(dpk.jac);
    
    err = dpGrad-dpGradNum;
    
    dpErrAbs1(k) = norm(   (err(mask))./gradScale(mask)  );
end
dpErrAbs1


dpErrAbs2 = zeros(nGrid,1);
for k = 1:nGrid
    [qok, qwk, qgk,pk] = deal(qoRandom(k), qwRandom(k), qgRandom(k), pRandom(k));
    
    dpGradNum =[
        (dpBeggsBrillJDJ(E, qok+qoP, qwk, qgk,pk)-dpBeggsBrillJDJ(E, qok-qoP, qwk, qgk,pk))/(2*qoP),...
        (dpBeggsBrillJDJ(E, qok, qwk+qwP, qgk,pk)-dpBeggsBrillJDJ(E, qok, qwk-qwP, qgk,pk))/(2*qwP),...
        (dpBeggsBrillJDJ(E, qok, qwk, qgk+qgP,pk)-dpBeggsBrillJDJ(E, qok, qwk, qgk-qgP,pk))/(2*qgP),...
        (dpBeggsBrillJDJ(E, qok, qwk, qgk,pk+ pP)-dpBeggsBrillJDJ(E, qok, qwk, qgk,pk- pP))/(2* pP),...
        ];
    
    
    [qok, qwk, qgk,pk] = initVariablesADI(qoRandom(k), qwRandom(k), qgRandom(k), pRandom(k));
    dpk = dpBeggsBrillJDJ(E, qok, qwk, qgk,pk);
    dpGrad = cell2mat(dpk.jac);
    
    err = dpGrad-dpGradNum;
    
    dpErrAbs2(k) = norm(   (err(mask))./gradScale(mask)  );
   % full([dpGrad',dpGradNum']);
end
dpErrAbs2





gScale = ((10*ft)^3/day);
gPert= gScale*1e-5;



nGrid = 1000;
qo = qoBounds(1) + (qoBounds(2)-qoBounds(1))*rand(1,1);
qw = qwBounds(1) + (qwBounds(2)-qwBounds(1))*rand(1,1);
p =   pBounds(1) + (pBounds(2) -pBounds(1)) *rand(1,1);

qgGrid =  linspace(qgBounds(1),qgBounds(2),nGrid);   



dp = zeros(nGrid,1);
dpErr = zeros(nGrid,1);
for k = 1:nGrid
    qg = qgGrid(k);
    
    dp(k) = dpBeggsBrillJDJ(E, qo, qw, qg,p);
    
    dpGradNum =(dpBeggsBrillJDJ(E, qo, qw, qg+gPert,p)-dpBeggsBrillJDJ(E, qo, qw, qg-gPert,p))/(2*gPert);
    
    
    [qok, qwk, qgk,pk] = initVariablesADI(qo, qw, qg, p);
    dpk = dpBeggsBrillJDJ(E, qok, qwk, qgk,pk);
    dpGrad = dpk.jac{3};
    
    
    dp(k) = dpk.val;
    
    
    dpErr(k) =norm(dpGradNum-dpGrad);
    
 
end


figure(1)
subplot(2,1,1)
plot(qgGrid/gScale,dp/barsa)
subplot(2,1,2)
plot(qgGrid/gScale,dpErr/barsa)

% 