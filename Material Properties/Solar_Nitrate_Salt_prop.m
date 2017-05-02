%This function will calculate the thermophysical properties of Solar Salt
%(NaNO3 ? KNO3)at a given temperature in Celsius and given pressure in bar.
%All correlations for Solar Salt obtained from 'Molten salts database for 
%energy applications' by R. Serrano-López,J Fradera, & S Cuesta-López 
%September 17, 2013. https://arxiv.org/pdf/1307.7343.pdf
function [mu_l,Cp_l,k_l,rho_l,Pr_l] = Solar_Nitrate_Salt_prop(T_l_avg) %Temp is in Celsius
T_l_avg_K=T_l_avg+273.15; %Converts temperature from Celsius to Kelvin
rho_l=2263.641-0.636*T_l_avg_K; %[kg/m^3]
mu_l=0.07543937-2.77*10^-4*T_l_avg_K+3.49*10^-7*T_l_avg_K^2-1.47*10^-10*T_l_avg_K^3; %[kg/m·s]
k_l=0.45; %[W/m·K]
Cp_l=1396.044+.172*T_l_avg_K;%[J/kg·K]
Pr_l=Cp_l*mu_l/k_l; %dimensionless