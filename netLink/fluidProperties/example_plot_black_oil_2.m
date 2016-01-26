% Script file example_plot_black_oil_2
%
% Creates plots of local (downhole) oil and gas rates as a function of pressure p at a fixed
% temperature T, for a given unit surface oil rate and the corresponding surface associated gas
% rate. The black oil properties are identical to those of example_plot_black_oil_1 and are
% computed with the aid of the Standing black oil correlations.     
%
% JDJ, 07-05-13, last revised 10-05-13
%

clear all
close all

% Set black oil and water data:
R_sb = 100; % gas_oil ratio at bubble point, m^3/m^3
rho_g_sc = 0.95; % gas density at standard conditions, kg/m^3
rho_o_sc = 850; % oil density at standard conditions, kg/m^3 
rho_w_sc = 1050; % water density at standard conditions, kg/m^3
oil = 1; % black oil; parameters computed with the aid of Standing correlations

% Set temperature:
T = 60; % temperature, deg. C

% Set surface rates:
q_o_sc = 1; % unit surface oil rate, m^3/s
q_g_sc = R_sb * q_o_sc; % corresponding surface associated gas rate, m^3/s
q_w_sc = 0; % zero surface water rate, m^3/s

% Create input vectors:
q_sc = [q_g_sc,q_o_sc,q_w_sc]; % flow rates at standard conditions, m^3/s
rho_sc = [rho_g_sc,rho_o_sc,rho_w_sc]; % densities at standard conditions, kg/m^3

% Create plot data:
n_step = 500; % number of steps in plot,-
Delta_p = 1e5; % pressure increment, Pa
p = 0; % initial pressure, Pa
results = zeros(n_step,3); % initializes matrix of plot data
for i = 1:n_step
    p = p + Delta_p; % pressure, Pa
    [q,rho] = local_q_and_rho(oil,p,q_sc,R_sb,rho_sc,T); % computes local rates and densities
    results(i,1) = p;
    results(i,2) = q(1); % local gas rate, m^3/s
    results(i,3) = q(2); % local oil rate, m^3/s
end

% Create plots:
subplot(2,2,1);
plot(results(:,1),results(:,3));
xlabel('Pressure\it p ,\rm Pa')
ylabel('Local oil rate\it q_o ,\rm m^3/s')
axis([0,5e7,1,1.4]) % sets lower and upper axis limits 
grid on
subplot(2,2,2);
plot(results(:,1),results(:,2));
xlabel('Pressure\it p ,\rm Pa')
ylabel('Local gas rate\it q_g ,\rm m^3/s')
axis([0,5e7,0,150])  
grid on
subplot(2,2,4);
semilogy(results(:,1),results(:,2));
xlabel('Pressure\it p ,\rm Pa')
ylabel('Local gas rate\it q_g ,\rm m^3/s')
axis([0,5e7,1e-3,1.5e2])  
grid on

