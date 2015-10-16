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

v1 = newVertex(1, -1,-1);    
%v1.pressure = 60*barsa;

v2  = newVertex(2, -1, -1); 
v2.pressure = 20*barsa;
%v2.pressure = 800*psia;

e1 = newEdge(1, v1, v2, -1);
e1.units = 0; % METRIC=0, FIELD = 1,
e1.pipeline = newPipeline('diam', 139.7/(10^3)*meter^3, 'len', 3535 , 'ang', 90*(pi/180), 'temp',  15.56*(274.15));
e1.stream = newStream();

e1.qoE = 700*meter^3/day;
e1.qwE = 150*meter^3/day;
e1.qgE = 0*meter^3/day; 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%  Edge e2  %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

v3 = newVertex(3, -1,-1);    

v4  = newVertex(4, -1, -1); 
v4.pressure = 30*barsa;

e2 = newEdge(1, v3, v4, -1);
e2.units = 1; % METRIC=0, FIELD = 1,
e2.pipeline = newPipeline('diam', 139.7/(10^3)*meter^3, 'len', 5535 , 'ang', 90*(pi/180), 'temp',  15.56*(274.15));
e2.stream = newStream();

e2.qoE = 700*meter^3/day;
e2.qwE = 200*meter^3/day;
e2.qgE = 0*meter^3/day; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%  ADI Objects  %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Network Graph: G = (V,E) %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
E = [e1; e2];
V = [v1; v2; v3; v4];

%E = [e1];
%V = [v1; v2];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Initing ADI Variables %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[E.qoE, E.qwE, E.qgE, V.pressure] = initVariablesADI(E.qoE, E.qwE, E.qgE, V.pressure);

%e1.pipeline = newPipeline('diam', 2.5*inch, 'len', 3535 , 'ang', 90*(pi/180), 'temp',  352.594); % temp = 175 F
%e1.stream = newStream('sg_gas', 41.6480048,'oil_dens', 799.321322 , 'water_dens', 0 , 'oil_visc', 2 ,'water_visc',0, 'gas_visc', 0.01331);
%e1.qoE = 0.00276019609;
%e1.qgE = 101940660;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Graph  0described as vectors of flows (qo, qw, qg) and pressures (p)  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Vout  = [];
for i=1:length(V)
    if V(i).id == 2 || V(i).id == 4
       Vout = [Vout; V(i)]; 
    end
end

[qo, qw, qg, p] = graph2FlowsAndPressures(Vout, E);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Calculating pressure drops %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dp = dpBeggsBrill(E, qo, qw, qg, p);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Validating pressure drops %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%[errorMax, errorJo, errorJw, errorJg] = testRunNetworkADI(E, qo, qw, qg, p);








