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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% Edge e1 %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

v1 = newVertex(1, -1,-1);   
v2  = newVertex(2, -1, -1); 
v2.pressure = 800*psia; % in Pa

e1 = newEdge(1, v1, v2, -1);
e1.units = 0; % METRIC=0, FIELD = 1,
e1.pipeline = newPipeline('diam', 2.5*inch, 'len', 500*meter , 'ang', degtorad(90), 'temp',  convtemp(175,'F','K'));
e1.stream = newStream();

e1.stream.gas_visc = 0.0131*(centi*poise); % viscosity in kilogram/(meter*second)
e1.stream.sg_gas = 0.65;
e1.stream.gas_dens = 2.6*pound/ft^3; % in kg/m^3

e1.stream.oil_visc = 2*(centi*poise); % viscosity in kilogram/(meter*second)
e1.stream.sg_oil = 0.35;
e1.stream.oil_dens =  49.9*pound/ft^3; % in kg/m^3

% e1.qoE = 317.97; % 2000 std
f = 800;
wcut = 0.4;

e1.qoE = -f*(meter^3/day);   % sm3/s
e1.qwE = -wcut*f*(meter^3/day);    % sm3/s

%e1.qgE = 46116*meter^3/day;  % GOR = 3393 scf/bbl = 604.3191 m3/m3
% e1.qgE = 28316.85; % in m3/d
% e1.qgE = 1*(10^6*ft^3/day)/(meter^3/day);
e1.qgE = 0*(meter^3/day);

E = [e1];
V = [v1; v2];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initing ADI Variables %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Graph  0described as vectors of flows (qo, qw, qg) and pressures (p)  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Vout  = [];
for i=1:length(V)
    if V(i).id == 2
       Vout = [Vout; V(i)]; 
    end
end

[qo, qw, qg, p] = graph2FlowsAndPressures(Vout, E);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Calculating pressure drops %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dp = simpleDp(E, qo, qw, qg, p);

dp = dpBeggsBrill(E, qo, qw, qg, p);
dp2 = simpleDp(E, qo, qw, qg, p);

dp/barsa
dp2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Validating pressure drops %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[errorMax, errorJo, errorJw, errorJg] = testRunNetworkADI(E, qo, qw, qg, p, 'qopert', 0.1*meter^3/day, 'qwpert', 0.1*meter^3/day, 'qgpert', 0);


