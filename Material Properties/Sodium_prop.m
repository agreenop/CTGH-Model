%This function will calculate the thermophysical properties of salt at a
%given temperature in Celsius and given pressure in bar.  It is based on Scarlat and Zweibaum, "Temperature-Dependent Thermophysical Properties
%for Fluoride Salts and Simulant Fluids", CIET-DESIGN-001-03, Rev. 03, 05-2-2013
function [mu_l,Cp_l,k_l,rho_l,Pr_l,nu_l] = Sodium_prop(Tl_avg) %Temp is in Celsius and Pressure in bar
Tl_avg=Tl_avg+273.15; %Convert to Kelvin
mu_l=813.9/(Tl_avg)-2.530; %kg/(m*s)
Cp_l=1732-1.266*Tl_avg+1.031*10^-3*Tl_avg^2-2.369*10^-7*Tl_avg^3; %J/(kg*K)
k_l=104+.0487*Tl_avg; % W/(m*K)
rho_l=1020-.247*Tl_avg; % kg/(m^3)
nu_l=mu_l/rho_l; %(m^2/s)
Pr_l=Cp_l*mu_l/k_l; %dimensionless