%This function will calculate the thermophysical properties of salt at a
%given temperature in Celsius and given pressure in bar.  It is based on Scarlat and Zweibaum, "Temperature-Dependent Thermophysical Properties
%for Fluoride Salts and Simulant Fluids", CIET-DESIGN-001-03, Rev. 03, 05-2-2013
function [mu_s,Cp_s,k_s,rho_s,Pr_s,nu_s] = Flibe_prop(Ts_avg) %Temp is in Celsius and Pressure in bar
mu_s=4.638*10^5/(Ts_avg)^2.79; %kg/(m*s)
Cp_s=2415.78; %J/(kg*K)
k_s=0.7662+0.0005*Ts_avg; % W/(m*K)
rho_s=2279.92-0.488*Ts_avg; % kg/(m^3)
nu_s=mu_s/rho_s; %(m^2/s)
Pr_s=Cp_s*mu_s/k_s; %dimensionless