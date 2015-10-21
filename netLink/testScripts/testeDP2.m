clear; clc;

% Required MRST modules
mrstModule add deckformat
mrstModule add ad-fi

addpath(genpath('../optimization/remso'));
addpath(genpath('../optimization/utils'));
addpath(genpath('../mrstDerivated'));
addpath(genpath('../netLink'));
addpath(genpath('../netLink/graphFunctions'));
addpath(genpath('../netLink/networkFunctions'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% Edge e1 %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pbar = 14.7*psia;

v1 = newVertex(1, -1,-1);   
v2  = newVertex(2, -1, -1); 
v2.pressure = 700; % bar

e1 = newEdge(1, v1, v2, -1);
e1.units = 0; % METRIC=0, FIELD = 1,
e1.pipeline = newPipeline('diam', 2.5*inch, 'len', 2737*ft , 'ang', degtorad(90), 'temp',  convtemp(175,'F','C'));
e1.stream = newStream();

e1.stream.gas_visc = 0.0131;
e1.stream.sg_gas = 0.65;
e1.stream.gas_dens = 2.6*pound/ft^3;

e1.stream.oil_visc = 2;
e1.stream.sg_oil = 0.35;
e1.stream.oil_dens =  49.9*pound/ft^3;

% e1.qoE = 317.97; % 2000 std
e1.qoE = -110;   %sm3/d
e1.qwE = -20;    % sm3/d

%e1.qgE = 46116*meter^3/day;  % GOR = 3393 scf/bbl = 604.3191 m3/m3
% e1.qgE = 28316.85; % in m3/d
e1.qgE = 0;

E = [e1];
V = [v1; v2];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initing ADI Variables %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[E.qoE, E.qwE, E.qgE, V.pressure] = initVariablesADI(E.qoE, E.qwE, E.qgE, V.pressure);

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
dp = dpBeggsBrill(E, qo, qw, qg, p);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Validating pressure drops %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[errorMax, errorJo, errorJw, errorJg] = testRunNetworkADI(E, qo, qw, qg, p, 'qgpert', 0);


