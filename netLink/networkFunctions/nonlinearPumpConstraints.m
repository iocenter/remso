function [ freq, qf, qpump_min, qpump_max, dhf, dpf] = nonlinearPumpConstraints( netSol, fref, numStages, qlMin, qlMax)
% nonlinearPumpConstraints: compute the nonlinear constraints of the pump
% operating envelope
%% compute local gas and liquid flor for pipe conditions
equipEdges = getEdge(netSol, netSol.Eeqp);

vin = vertcat(equipEdges.vin);
vout = vertcat(equipEdges.vout);
    
dpf = netSol.pV(vin)-netSol.pV(vout); % dp in the pumps

[~, ~, ~,~, oil, rho_sc, ~, ~ , ~, T_in,~] = wrapperJDJ(equipEdges);
q_g_sc = netSol.qg(vertcat(equipEdges.id));
q_o_sc = netSol.qo(vertcat(equipEdges.id));
q_w_sc = netSol.qw(vertcat(equipEdges.id));
pInt = netSol.pV(vertcat(equipEdges.vin)); % pump intake pressure
T = T_in; % pump intake temperature
hasSurfaceGas = false;
Z = [];
R_sb = zeros(size(double(q_o_sc)));

[q_g,q_o,q_w,~,rho_o,rho_w,~] = local_q_and_rho(oil,pInt,q_g_sc,q_o_sc,q_w_sc,R_sb,rho_sc,T,hasSurfaceGas,Z);
if any(q_g >= 1e-03*meter^3/day)
    warning('Gas flow appeared at pump local conditions.');
end
qf = q_o + q_w;
mixtureDen = (q_o.*rho_o + q_w.*rho_w)./qf;

dhf= pump_dh(dpf, mixtureDen); % dh in the pumps

%% TODO: include parameter with the reference frequency
freq = pump_eq_system_explicit(qf, dhf, fref, numStages);  % solves a system of equations to obtain frequency, flow and dh at 60Hz

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Flow Constraint in Equipment  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
qpump_min = pump_rate(freq, qlMin, 60);
qpump_max = pump_rate(freq, qlMax, 60);

end

