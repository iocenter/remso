function dp =  dpBeggsBrill(E, qoE, qwE, qgE, pV, hasSurfaceGas, f_work, Z)
%% dp: calculates the pressure drop for the average pressure pV
%%     the functions is implemented for SI units (input).

if nargin < 7
    f_work = [];
    Z = [];
end
if nargin < 6
    hasSurfaceGas = true;
end


%% Validate the function
% Compare against simpleDP and PIPESIM for a given pipeline.

%     units = struct('METRIC',0, 'FIELD', 1);

nE = length(E);
str = vertcat(E.stream);

%% Producer flows are negative, and injection flows are positive.
emptyADI = getEmptyADI(qoE, qwE, qgE, pV);
dp = emptyADI;
zfactor = emptyADI;

total_rate = qgE + qoE + qwE;
flag_rate = abs(double(total_rate)) >= 1e-06*meter^3/day;
flag_gas = abs(double(qgE)) > 0;  %%TODO: determine what is a very little gas rate. Is it correct to set it to zero ?

if any(~flag_rate)
    if any(flag_rate)
        dp(flag_rate) =  dpBeggsBrill(E(flag_rate), qoE(flag_rate), qwE(flag_rate), qgE(flag_rate), pV(flag_rate), hasSurfaceGas, f_work(flag_rate), Z(flag_rate));
    end
    return
end

pipes = vertcat(E.pipeline);
diameters = vertcat(pipes.diam);
temperatures = vertcat(pipes.temp);
angles = vertcat(pipes.ang);

grav = norm(gravity);                                       % gravitational constant (m/s^2)
if any(flag_gas)
    strGas = vertcat(str.sg_gas);
    
    zfactor(flag_gas) = zFactor(strGas(flag_gas), temperatures(flag_gas), pV(flag_gas));     % gas z-factor
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Flow Velocities %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
vsl         = superficialLiquidVelocity(qoE, qwE,diameters);           % superficial liquid velocity returns in m/s
vsg         = emptyADI;
if any(flag_gas)
    vsg(flag_gas)   = superficialGasVelocity(qgE(flag_gas), pV(flag_gas), zfactor(flag_gas), diameters(flag_gas), temperatures(flag_gas)); % superficial two phase velocity in m/s
end
vm         = vsl  + vsg;

froude_num = vm.^2./diameters./grav;                                 % froude number
liquid_content = vsl./vm;

%% if the liquid content is 0, changing it to something small
%     liquid_content = max(liquid_content,  1e-8.*ones(nE,1));
if any(liquid_content < 1e-8.*ones(nE,1))
    warning('Liquid flow is approaching to zero.');
end

%% calculation of the correction factor
%% TODO: check if pressure should affect the density
den_l = liquidDensity(qoE, qwE, str);                          %% liquid density in SI
den_g = emptyADI;
if  any(flag_gas)
    den_g(flag_gas) = gasDensity(gasSpecificGravity(str(flag_gas)), temperatures(flag_gas), pV(flag_gas), zfactor(flag_gas));    %% gas density in kg/m^3
end


holdup    = emptyADI + ones(nE,1);
if any(flag_gas)
    holdup(flag_gas) = calcHoldup(liquid_content(flag_gas),froude_num(flag_gas),den_l(flag_gas),den_g(flag_gas),vsl(flag_gas),emptyADI(flag_gas),grav,angles(flag_gas));
end

%% two phase density
den_s = den_l .* holdup + den_g .* (1 - holdup);

%% calculating pressure drop due to elevation change
dp_el = den_s.*grav.*sin(angles);    %% pressure gradient due to elevation change

%% calculating friction factor
vis_g = emptyADI;
if any(flag_gas)
    vis_g(flag_gas) = gasViscosity(str(flag_gas), temperatures(flag_gas), pV(flag_gas), zfactor(flag_gas));
end

den_ns = den_l .* liquid_content + den_g .* (1 - liquid_content);      %% no-slip density
vis_ns = liquidViscosity(qoE, qwE, str) .* liquid_content + vis_g .* (1 - liquid_content);       %% no-slip viscosity

re_ns = (den_ns .* abs(vm) .* diameters)./vis_ns;     %% no-slip reynolds number


%% reynolds_threshold = 10^(3.8215/4.5223) ~≃ 7
reynolds_threshold = 7;
condR = (re_ns > reynolds_threshold);
fn = re_ns;
if any(condR)
    fn(condR) = 1./(2 .* log((re_ns(condR) ./ (4.5223 .* log(re_ns(condR))./log(10) - 3.8215)))./log(10)).^2;    %% no-slip friction factor
end
if any(~condR)
    fn(~condR) = 0.0056 + 0.5./(re_ns(~condR)).^(0.32); % simple calculation for fn
end
%% TODO: consider including the friction factor for laminar flow... see page 6 of the book James Brill multiphase flow
%% TODO: check the effect of roughness


cond_fn0 = fn < 0;
if any(cond_fn0)
    warning('Friction factor is negative. Setting up to zero.')
    fn(cond_fn0) = 0;
end
y = liquid_content ./(holdup.^2);

cond_y = (y > 1.0 & y < 1.2);
s_term = y; % initialization, does not require to be y.
if any(cond_y)
    s_term(cond_y) = log(2.2.*y(cond_y) -1.2);
end
if any(~cond_y)
    s_term(~cond_y) = log(y(~cond_y)) ./ (-0.0523 + 3.182 .* log(y(~cond_y)) - 0.8725 .* (log(y(~cond_y)).^2) + 0.01853 .* (log(y(~cond_y)).^4));
end

ftp = fn .* exp(s_term);      %% the friction factor

%% calculating pressure drop due to friction
dp_f = -(ftp.*den_ns.*abs(vm).*vm)./(2.*diameters); %% James P. Brill, H. Dale Beggs Two-Phase  Flow in Pipes


%% calculating acceleration term
ek = emptyADI; %% it has impact only when there is gas flowing in the pipe. if there is no gas then vsg==0
if any(flag_gas)
    ek(flag_gas) = (den_s(flag_gas) .* vm(flag_gas) .* vsg(flag_gas)) ./ pV(flag_gas);
end

%% calculating total pressure drop gradient

%% pressure gradient
dp = (dp_f + dp_el) ./ (1 - ek);  %% total pressure drop per length of pipe (in Pa./ m)
end

function [holdup] = liqholdupTransitionFlow(liquid_content, nlv, froude_num, angles, l2, l3, payne_cor)
frac = (l3-froude_num)./(l3-l2);

%% SEGREGATED flow constants
a_seg = 0.98;
b_seg = 0.4846;
c_seg = 0.0868;

%% INTERMITTENT flow constants
a_int = 0.845;
b_int = 0.5351;
c_int = 0.0173;

%% horizontal holdups
hz_holdup_seg = (a_seg.*(liquid_content.^b_seg))./(froude_num.^c_seg);
hz_holdup_int = (a_int.*(liquid_content.^b_int))./(froude_num.^c_int);

%% horizontal holdup is the liquid content if smaller.
%% -(max(-x,-y)) = min(x,y) for positive integers
hz_holdup_seg = max(hz_holdup_seg, liquid_content);
hz_holdup_int = max(hz_holdup_int, liquid_content);

%% correction factors
cor_seg = liquid_content.*0;
cor_int = liquid_content.*0;

cond_ang = (angles < 0);
if any(cond_ang)
    d_ang = 4.7;
    e_ang = -0.3692;
    f_ang =  0.1244;
    g_ang = -0.5056;
    
    cor_seg(cond_ang) = (1-liquid_content(cond_ang)).*log(d_ang.*(liquid_content(cond_ang).^e_ang).*(nlv(cond_ang).^f_ang).*(froude_num(cond_ang).^g_ang));
    cor_int(cond_ang) = cor_seg(cond_ang);
    
end

if any(~cond_ang)
    d_seg = 0.011;
    e_seg = -3.768;
    f_seg = 3.539;
    g_seg = -1.614;
    
    d_int = 2.96;
    e_int = 0.305;
    f_int = -0.4473;
    g_int = 0.0978;
    
    cor_seg(~cond_ang) = (1 - liquid_content(~cond_ang)) .* log(d_seg.*(liquid_content(~cond_ang).^e_seg).*(nlv(~cond_ang).^f_seg).*(froude_num(~cond_ang).^g_seg));
    cor_int(~cond_ang) = (1 - liquid_content(~cond_ang)) .* log(d_int.*(liquid_content(~cond_ang).^e_int).*(nlv(~cond_ang).^f_int).*(froude_num(~cond_ang).^g_int));
end
cor_seg = max(cor_seg, cor_seg*0);
cor_int = max(cor_int, cor_int*0);

psi_seg = 1 + cor_seg .* (sin(1.8.*angles) - 0.333 .* (sin(1.8.*angles)).^3 );
psi_int = 1 + cor_int .* (sin(1.8.*angles) - 0.333 .* (sin(1.8.*angles)).^3);

holdup = payne_cor .* (frac .* (hz_holdup_seg .* psi_seg) + (1 - frac) .* (hz_holdup_int .* psi_int));
%        holdup = (frac .* (hz_holdup_seg .* psi_seg) + (1 - frac) .* (hz_holdup_int .* psi_int));
end

%% gasSpecificGravity: returns the gas specific gravity
function sg = gasSpecificGravity(str)
sg = vertcat(str.sg_gas);
end

%% superficialGasVelocity: calculates the superficial gas velocity (vsg)
function gv = superficialGasVelocity(qgE, pres, zfac, diam, temp)
gas_rate_surface = qgE;         % gas rate in sm3/s
% surface conditions
Tsc = convtemp(60,'F','K');     % surface temperature in K
Psc = 14.7*psia;                % surface pressure in Pa

pres = pres./barsa;
Psc  = Psc./barsa;

% pipe conditions
A = pi.*((diam./2)).^2;     % pitpe area in m^2

num = gas_rate_surface.*zfac.*temp.*Psc;
den = A.*Tsc.*pres;

gv = num./den;

end


%% superficialLiquidVelocity: calculates the superficial liquid velocity (vsl)
function liqVeloc =  superficialLiquidVelocity(qo, qw, diam)

liquid_rate = (qo + qw);              % liquid rate in sm3/s

area =  pi.*((diam./2)).^2;                             % pipe radius in m^2

liqVeloc = liquid_rate./area;                 % in m/s
end

%% liqDens: calculates the liquid density (oil and water)
function liqDens = liquidDensity(qoE, qwE, str)
oil_rate = qoE;
water_rate = qwE;

if any(abs(oil_rate + water_rate) < 1.e-10*meter^3/day)
    warning('Liquid rate approaching to zero. Impossible to calculate liquid density.');
end
oil_mdensity = oil_rate.*vertcat(str.oil_dens);
water_mdensity =  water_rate.*vertcat(str.water_dens);


liqDens = (oil_mdensity + water_mdensity)./(oil_rate + water_rate);   % liquid density in kg/m^3
end
function zfac = zFactor(sg, t, p)
%% gasZFactor: calculates gas z-factor
% t - temperature
% p - pressure
% sg - gas specific gravity
zfac = Z_factor_DAK_direct(p,sg,t);
end

%% Calculates the gas density at pipe conditions
function denGas = gasDensity(sg, temp, pres, zfac)
air_molecular_weight = 28.97*gram/(kilogram);                                    % in kg/mol
ideal_gas_cons = 8.3144598;                                                      % J/mol K

Mg = sg.*air_molecular_weight;                                                   % molecular weight of gas

denGas  = (pres .* Mg)./(ideal_gas_cons.*zfac.*temp);                    % in kg/m.^3
end

%% Calculates the gas-liquid surface tension
function db = surfaceTension(gas_density, liquid_density)
%% converting to original unity in resopt to perform calculations
%% Reference for this correlation can be found in a report written by prof. Curtis H. Whitson in November 20, 1990 entitled 'Minimum Lift Calculation for Gas Wells'.
%% The correlation was originally proposed by Ramey (SPE 4429).

gas_density_field = gas_density./(pound/(ft^3));               % gas density in lb/ft^3
liq_density_field = liquid_density./(pound/(ft^3));            % liq density in lb/ft^3


db_field = liq_density_field - gas_density_field;
db_field = 15.0 + 0.91.*db_field;                              % surface tension in pound/s^2

db_field = db_field*(pound/second^2);                                % in dyne/cm

db = db_field.*(dyne)./(centi*meter);                          % surface tension in J / m
end


%% gasViscosity: calculates the gas viscosity
function viscGas = gasViscosity(str, t, p, z)
% E - pipeline
% t - temperatures
% p - pressure
% z - zfactor
% Using correlation found in the phd thesis of Aleks Juell.
% The correlation have first appeared in the SPE paper
% 'The Viscosity of Natural Gases' by Anthon L. Lee, Mario H. Gonzale
% In the experimental settings, the density is in g/cm3,
% the temperature is in R and the pressure in psia.


t_r = convtemp(t, 'K', 'R'); % temperature in R
sg = gasSpecificGravity(str);


den_gas_SI = gasDensity(sg, t, p, z); %% gas density in kg/m^3
den_gas = den_gas_SI/(gram/(centi*meter)^3);  % in g/cm^3

%     calculation to obtain density from PV = NRZT formulae
%     den_gas = vertcat(str.gas_dens)/(gram/(centi*meter)^3);

air_molecular_weight = 28.97;                                   % in /mol
Mg = sg.*air_molecular_weight;                                  % molecular weight of gas


A1 = ((9.379 + 0.01607.*Mg).*t_r.^1.5)./(209.2  + 19.26.*Mg + t_r);
A2 = (3.448 + 986.4./t_r + 0.01009.*Mg);
A3 = 2.447 - 0.2224.*A2;

viscGas = A1.*exp(A2.*(den_gas.^A3));  % in micropoise

viscGas = viscGas*(micro*poise); %% kilogram per meter-second (SI)
end


%% liquidViscosity: calculates the liquid viscosity
function liqVisc = liquidViscosity(qo, qw, str)
oil_rate = qo;
water_rate = qw;

cond_lowRate = abs(oil_rate + water_rate) < 1e-8*meter^3/day;
if any(cond_lowRate)
    warning('Liquid rate approaching to zero. Impossible to calculate liquid viscosity.');
    oil_rate(cond_lowRate) = 1.0;
end

liqVisc = ( oil_rate.*vertcat(str.oil_visc) +water_rate.*vertcat(str.water_visc) )./(oil_rate + water_rate);
end

function emptyADI = getEmptyADI(qoE, qwE, qgE, pV)

if isa(qoE,'ADI')
    emptyADI = qoE*0;
elseif isa(qwE,'ADI')
    emptyADI = qwE*0;
elseif isa(qgE,'ADI')
    emptyADI = qgE*0;
elseif isa(pV,'ADI')
    emptyADI = pV*0;
else
    emptyADI = zeros(size(pV));  %% it is a double
end

end


function holdup = calcHoldup(liquid_content,froude_num,den_l,den_g,vsl,emptyADI,grav,angles)

flow_regime = struct('SEGREGATED', 0, 'TRANSITION', 1, 'INTERMITTENT', 2, 'DISTRIBUTED', 3, 'UNDEFINED', 4);

nE = numel(double(liquid_content));

%% constants used in formula for evaluating the flow regime
l1 = 316.*liquid_content.^(0.302);
l2 = 0.0009252.*liquid_content.^(-2.4684);
l3 = 0.10.*liquid_content.^(-1.4516);
l4 = (0.5).*(liquid_content.^(-6.738));

%% conditional vectors which are used to determine flow regimes
cond_seg = (liquid_content < 0.01 & froude_num < l1) | (liquid_content>=0.01 & froude_num < l2);
cond_trans = (liquid_content >= 0.01 & froude_num >= l2 & froude_num <= l3);
cond_dist = (liquid_content < 0.4 & froude_num >= l1) | (liquid_content >= 0.4 & froude_num > l4);
cond_interm = (liquid_content >= 0.01 & liquid_content < 0.4 & froude_num > l3 & froude_num <= l1) | ...
    (liquid_content >= 0.4 & froude_num > l3 & froude_num <= l4) | ...
    (~cond_seg & ~cond_trans & ~cond_dist);  % assuming intermittent flow when it is not possible to determine.  TODO: Is this cheating?


surface_tens = surfaceTension(den_g, den_l);          %% gas - liquid surface tension in SI   %%TODO: confirm the correctness of the correlation when there is no gas flowing

nlv = vsl .* (den_l./(grav.*surface_tens)).^(0.25);         %% liquid velocity number

%% Payne correction factor to holdup  (Page 43 of the Book 'James Brill, Multiphase Flow in Wells, No 17')
payne_cor = 0.924.*ones(nE,1);

% Beggs % Brill holdup constants
regime = zeros(nE,1);
cor    = emptyADI;
psi    = emptyADI + ones(nE,1);

a = zeros(nE,1);
b = zeros(nE,1);
c = zeros(nE,1);

d = zeros(nE,1);
e = zeros(nE,1);
f = zeros(nE,1);
g = zeros(nE,1);

cond_ang = (angles < zeros(nE,1));
if any(cond_ang) % Negative Angle
    d(cond_ang) = 4.7;
    e(cond_ang) = -0.3692;
    f(cond_ang) = 0.1244;
    g(cond_ang) = -0.5056;
    
    payne_cor(cond_ang) = 0.685;
end


% observe that d,e,f,g are reset when cond_ang is true
% why don't you follow the name conventions in the book 'James Brill, Multiphase Flow in Wells, No 17'
if any(cond_seg) % SEGREGATED Flow
    a(cond_seg) = 0.98;
    b(cond_seg) = 0.4846;
    c(cond_seg) = 0.0868;
    
    d(cond_seg & ~cond_ang) = 0.011;
    e(cond_seg & ~cond_ang) = -3.768;
    f(cond_seg & ~cond_ang) = 3.539;
    g(cond_seg & ~cond_ang) = -1.614;
    
    regime(cond_seg) =  flow_regime.SEGREGATED;
end

if any(cond_interm) % INTERMITTENT flow
    
    a(cond_interm) = 0.845;
    b(cond_interm) = 0.5351;
    c(cond_interm) = 0.0173;
    
    d(cond_interm & ~cond_ang) = 2.96;
    e(cond_interm & ~cond_ang) = 0.305;
    f(cond_interm & ~cond_ang) = -0.4473;
    g(cond_interm & ~cond_ang) = 0.0978;
    
    regime(cond_interm) = flow_regime.INTERMITTENT;
    
end
if any(cond_dist) % DISTRIBUTED flow
    
    a(cond_dist) = 1.065;
    b(cond_dist) = 0.5824;
    c(cond_dist) = 0.0609;
    
    d(cond_dist & ~cond_ang) = 0;
    e(cond_dist & ~cond_ang) = 0;
    f(cond_dist & ~cond_ang) = 0;
    g(cond_dist & ~cond_ang) = 0;
    
    regime(cond_dist) = flow_regime.DISTRIBUTED;
    
end
if any(cond_trans) % TRANSITION flow
    regime(cond_trans) = flow_regime.TRANSITION;
end
assert(all(or(cond_interm,or(cond_seg,or(cond_trans,cond_dist)))),'Couldn''t find flow regime');


cond_not_uphillAndDistorTrans = ~(((~cond_ang) & cond_dist) | cond_trans );

if any(cond_not_uphillAndDistorTrans)
    cor(cond_not_uphillAndDistorTrans) = (1 - liquid_content(cond_not_uphillAndDistorTrans)).*...
        log(d(cond_not_uphillAndDistorTrans).*...
        (liquid_content(cond_not_uphillAndDistorTrans).^e(cond_not_uphillAndDistorTrans)).*...
        (nlv(cond_not_uphillAndDistorTrans).^f(cond_not_uphillAndDistorTrans)).*...
        (froude_num(cond_not_uphillAndDistorTrans).^g(cond_not_uphillAndDistorTrans)));
end
%if any(~cond_not_uphillAndDistorTrans)  it is already defined as zero!
%     cor(~cond_not_uphillAndDist) = 0
%end

%% Horizontal liquid holdup, hl(0), except for TRANSITION flow,
%% which will be later calculated as an interpolation of SEGREGATED and
%% INTERMITTENT flows.
holdup = emptyADI;
if any(~cond_trans)
    hz_holdup = (a(~cond_trans).*(liquid_content(~cond_trans).^b(~cond_trans)))./(froude_num(~cond_trans).^c(~cond_trans));
    
    
    %% the horizontal holdup is equal to liquid content if smaller
    hz_holdup = max(hz_holdup, liquid_content(~cond_trans));
    
    %% checking if correction >= 0
    cor(~cond_trans) = max(cor(~cond_trans), emptyADI(~cond_trans));
    %   condCor = (cor < 0);
    %   if any(condCor)
    %        disp '### Warning ###'
    %        disp 'From: Beggs & Brill 1973'
    %        disp 'The calculated correction factor, C, is negative... '
    %        disp 'Resetting to 0.0'
    %        warning('The calculated correction factor, C, is negative... Resetting to 0.0');
    %       cor(condCor) = zeros(sum(condCor),1);
    %   end
    
    %% liquid holdup correction for inclination
    if any(cond_not_uphillAndDistorTrans)
        psi(cond_not_uphillAndDistorTrans) = ones(length(angles(cond_not_uphillAndDistorTrans)),1) + cor(cond_not_uphillAndDistorTrans).*(sin(1.8.*angles(cond_not_uphillAndDistorTrans)) - 0.333.*(sin(1.8.*angles(cond_not_uphillAndDistorTrans)).^3));
    end
    % if not psi is already 1
    
    %% if otherwise psi is equal to 1 already
    
    %    holdup = hz_holdup.*psi; % liquid holdup corrected for inclination
    holdup(~cond_trans) = payne_cor(~cond_trans).*hz_holdup.*psi(~cond_trans); % liquid holdup corrected for inclination
    
    %% if trasition regime, the liquid holdup is a mix of segregated and intermmitent
end
if any(cond_trans)
    holdup(cond_trans) = liqholdupTransitionFlow(liquid_content(cond_trans), nlv(cond_trans), froude_num(cond_trans), angles(cond_trans), l2(cond_trans), l3(cond_trans), payne_cor(cond_trans));
end

%% corrects the holdup if wrong
holdup = -(max(-holdup, -ones(length(holdup),1)));


end
