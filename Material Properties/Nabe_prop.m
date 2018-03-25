%This function will calculate the thermophysical properties of nabe salt at a
%given temperature in Celsius and given pressure in bar.
function [mu_l,Cp_l,k_l,rho_l,Pr_l] = Nabe_prop(T_l_avg) %Temp is in Celsius
T_l_avg_K=T_l_avg+273.15; %Convert temperature to Kelvin
rho_l=2270-0.37*T_l_avg; % kg/(m^3)
mu_l=0.0346*10^-3*exp(5165/T_l_avg_K); %kg/(m*s)
Cp_l=2175.68; %J/(kg*K)
k_l=0.5222+0.0005*T_l_avg; % W/(m*K)
Pr_l=Cp_l*mu_l/k_l; %dimensionless