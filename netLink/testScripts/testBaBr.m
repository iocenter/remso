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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% Edge e1 %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

v1 = newVertex(1, -1,-1);   
v2  = newVertex(2, -1, -1); 
% v2.pressure = 800*psia; % in Pa

% v2.pressure = 720*psia; % in Pa

v2.pressure = 300*barsa;

e1 = newEdge(1, v1, v2, -1);
e1.units = 0; % METRIC=0, FIELD = 1,
% e1.pipeline = newPipeline('diam', 2.5*inch, 'len', 500*meter , 'ang', degtorad(90), 'temp',  convtemp(175,'F','K'));

e1.pipeline = newPipeline('diam', 0.249*ft, 'len', 1*ft , 'ang', degtorad(90), 'temp',  convtemp(175,'F','K'));

e1.stream = newStream();
% e1.stream.gas_visc = 0.0131*(centi*poise); % viscosity in kilogram/(meter*second)
e1.stream.gas_visc = 0.018*(centi*poise); % viscosity in kilogram/(meter*second)
e1.stream.sg_gas = 0.65;
% e1.stream.gas_dens = 2.6*pound/ft^3; % in kg/m^3
e1.stream.gas_dens = 2.84*pound/ft^3; % in kg/m^3

% e1.stream.oil_visc = 2*(centi*poise); % viscosity in kilogram/(meter*second)
e1.stream.oil_visc = 18*(centi*poise); % viscosity in kilogram/(meter*second)
e1.stream.sg_oil = 0.35;
% e1.stream.oil_dens =  49.9*pound/ft^3; % in kg/m^3
e1.stream.oil_dens =  56.6*pound/ft^3; % in kg/m^3

f = 319.68;
wcut = 0.15;

e1.qoE = -(1-wcut)*f*(meter^3/day);   % sm3/s
e1.qwE = -wcut*f*(meter^3/day);       % sm3/s
e1.qgE = -2583.84*(meter^3/day);             % sm3/s


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
e2.qgE = -10^6*(ft^3/day);             % sm3/s

E = [e1; e2];
V = [v1; v2; v3; v4];

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
dp = dpBeggsBrill(E, qo, qw, qg, p);
dp/barsa;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Validating pressure drops %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[errorMax, errorJo, errorJw, errorJg, errorJp] = testRunNetworkADI(E, qo, qw, qg, p, 'qopert', 1e-06*meter^3/day, 'qwpert', 1e-06*meter^3/day, 'qgpert',  1-06*meter^3/day, 'poutPert', 1e-03*barsa);


