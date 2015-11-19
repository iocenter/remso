% Script file example_plot_black_oil_1
%
% Creates plots of black oil parameters as a function of pressure p at a fixed temperature T
% using the Standing black oil correlations. 
%
% JDJ, 15-02-05, last revised 10-05-13
%

clear all
close all

% Set black oil data:
R_sb = 100; % gas_oil ratio at bubble point, m^3/m^3
rho_g_sc = 0.95; % gas density at standard conditions, kg/m^3
rho_o_sc = 850; % oil density at standard conditions, kg/m^3

% Set temperature:
T = 60; % temperature, deg. C

% Create plot data:
n_step = 500; % number of steps in plot,-
Delta_p = 1e5; % pressure increment, Pa
p = 0; % initial pressure, Pa
results = zeros(n_step,4); % initializes matrix of plot data
for i = 1:n_step
    p = p + Delta_p; % pressure, Pa
    [B_g,B_o,R_s] = black_oil_Standing(p,R_sb,rho_g_sc,rho_o_sc,T); % computes black oil param.
    results(i,1) = p;
    results(i,2) = B_g;
    results(i,3) = B_o;
    results(i,4) = R_s;
end

% Create plots:
subplot(2,2,1);
plot(results(:,1),results(:,3));
xlabel('Pressure\it p ,\rm Pa')
ylabel('Oil formation volume factor\it  B_o ,\rm m^3/m^3')
axis([0,5e7,1,1.4]) % sets lower and upper axis limits 
grid on
subplot(2,2,2);
plot(results(:,1),results(:,2));
xlabel('Pressure\it p ,\rm Pa')
ylabel('Gas formation volume factor\it  B_g ,\rm m^3/m^3')
axis([0,5e7,0,1]) 
grid on
subplot(2,2,3);
plot(results(:,1),results(:,4));
xlabel('Pressure\it p ,\rm Pa')
ylabel('Solution gas-oil ratio\it R_s ,\rm m^3/m^3')
axis([0,5e7,0,120]) 
grid on
subplot(2,2,4);
semilogy(results(:,1),results(:,2));
xlabel('Pressure\it p ,\rm Pa')
ylabel('Gas formation volume factor\it  B_g ,\rm m^3/m^3')
axis([0,5e7,1e-3,1e0]) 
grid on

