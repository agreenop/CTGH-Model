%This function will calculate and interpolate the thermophysical properties
%of water based on the data for saturated water in F. Incropera and D. De Witt,
%"Fundamentals of Heat Transfer," 5th Ed., Wiley, 1990. Table A.6, pg. 860
%based on the temperature given in Celsius.  It is assumed that
%thermophysical properties of saturated water are approximately the same for 
%non-saturated liquid water.
function [mu_l,Cp_l,k_l,rho_l,nu_l,Pr_l] = Water_prop(T_l_avg) %Temp is in Celsius
T_l_avg_K=T_l_avg+273.15;
water=csvread('water_data.csv',1,0); %Water data from Incropera
rho_l=interp1(water(:,1),water(:,2),T_l_avg_K); %[kg/m^3]
Cp_l=1000*interp1(water(:,1),water(:,3),T_l_avg_K); %[J/kg*K]
mu_l=interp1(water(:,1),water(:,4),T_l_avg_K); %[N*s/m^2]
k_l=interp1(water(:,1),water(:,5),T_l_avg_K); %[W/m*K]
nu_l=mu_l/rho_l;%(m^2/s)
Pr_l=Cp_l*mu_l/k_l; %dimensionless