function p_b = pres_bub_volatile_oil(T)
% function p_b = pres_bub_volatile_oil(T)
%
% Determines the bubble point pressure p_b of a tabulated volatile oil at a given temperature. 
%
% p_b = bubble point pressure, Pa
% T = temperature, deg. C
%
% JDJ, 08-05-13, last revised 27-10-13

% Load volatile oil data:
file_name = 'vol_oil_table_01';
load(file_name,'vol_oil');

% Step through temperature values to find nearest temperature:
i = 1;
T_lo = vol_oil(1,1,1); % lowest T value in table, deg. C
T_tab = T_lo;
while T_tab < T
    i = i + 1;
    T_tab = vol_oil(i,1,1);
end

% Step through pressure values until no further increase in R_s occurs:
j = 1;
R_s_lo = vol_oil(1,1,6); % lowest R_s value in table, m^3/m^3
R_s_tab = R_s_lo;
R_s_old = R_s_lo - eps; % to get loop started
while R_s_tab > R_s_old
    j = j + 1;
    R_s_old = R_s_tab;
    R_s_tab = vol_oil(i,j,6);
end
p_b = vol_oil(i,j,2);
