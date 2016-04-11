function f = gas_visc_ratio_Dempsey(p_pr,T_pr)
% function f = gas_visc_ratio_Dempsey(p_pr,T_pr)
%
% Calculates the ratio f between the gas viscosity at any pressure and the viscosity at
% atmosperic pressure for a given pseudo-reduced pressure and temperature. 
%
% Use is made of an expression of Dempsey (1965) to approximate the correlation of Carr,
% Kobayashi and Burrows (1954). 
%
% f = gas viscosity ratio = mu_g / mu_g_p_sc, -
% p_pr = pseudo-reduced pressure, -
% T_pr = pseudo-reduced temperature, -
%
% JDJ, 04-10-02, last revised 10-05-13
%{
%Changes by Thiago and Codas
%Make the function compatible with ADI objects
%Vectorization
%}

 a0 = -2.46211820e-00;
 a1 =  2.97054714e-00;
 a2 = -2.86264054e-01;
 a3 =  8.05420522e-03;
 a4 =  2.80860949e-00;
 a5 = -3.49803305e-00;
 a6 =  3.60373020e-01;
 a7 = -1.04432413e-02;
 a8 = -7.93385684e-01;
 a9 =  1.39643306e-00;
a10 = -1.49144925e-01;
a11 =  4.41015512e-03;
a12 =  8.39387178e-02;
a13 = -1.86408848e-01;
a14 =  2.03367881e-02;
a15 = -6.09579263e-04;

help01 =            a0 +  a1*p_pr +  a2*p_pr.^2 +  a3*p_pr.^3 ;
help02 = T_pr   .* ( a4 +  a5*p_pr +  a6*p_pr.^2 +  a7*p_pr.^3);
help03 = T_pr.^2 .* ( a8 +  a9*p_pr + a10*p_pr.^2 + a11*p_pr.^3);
help04 = T_pr.^3 .* (a12 + a13*p_pr + a14*p_pr.^2 + a15*p_pr.^3);

f = exp(help01+help02+help03+help04) ./ T_pr;
