%This function will calculate the thermophysical properties of salt at a
%given temperature in Celsius and given pressure in bar.  It is based on Scarlat and Zweibaum, "Temperature-Dependent Thermophysical Properties
%for Fluoride Salts and Simulant Fluids", CIET-DESIGN-001-03, Rev. 03, 05-2-2013
function [mu_l,Cp_l,k_l,rho_l,Pr_l] = Flibe_prop(T_l_avg) %Temp is in Celsius and Pressure in bar
mu_l=4.638*10^5/(T_l_avg)^2.79; %kg/(m*s)
Cp_l=2415.78; %J/(kg*K)
k_l=0.7662+0.0005*T_l_avg; % W/(m*K)
rho_l=2279.92-0.488*T_l_avg; % kg/(m^3)
Pr_l=Cp_l*mu_l/k_l; %dimensionless