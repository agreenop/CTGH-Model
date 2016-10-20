%This function will calculate the thermophysical properties of liquid sodium at a
%given temperature in Celsius.  These equations are based on "Database of thermophysical
%properties of liquid metal coolants for GEN-IV."
%www.iaea.org/inis/collection/NCLCollectionStore/_Public/43/095/43095088.pdf
function [mu_l,Cp_l,k_l,rho_l,Pr_l,nu_l] = Sodium_prop(T_l_avg) %Temp is in Celsius
T_l_avg=T_l_avg+273.15; %Convert to Kelvin
mu_l=813.9/(T_l_avg)-2.530; %kg/(m*s)
Cp_l=1732-1.266*T_l_avg+1.031*10^-3*T_l_avg^2-2.369*10^-7*T_l_avg^3; %J/(kg*K)
k_l=104+.0487*T_l_avg; % W/(m*K)
rho_l=1020-.247*T_l_avg; % kg/(m^3)
nu_l=mu_l/rho_l; %(m^2/s)
Pr_l=Cp_l*mu_l/k_l; %dimensionless