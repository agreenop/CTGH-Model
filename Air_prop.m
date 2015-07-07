%This function will calculate and interpolate the thermophysical properties
%of the air based on the data for air in F. Incropera and D. De Witt,
%"Fundamentals of Heat Transfer," 5th Ed., Wiley, 1990. Table A.4, pg. 852
%based on the temperature given in Celsius and pressure in bar.
function [rho_a,cp_a,mu_a,k_a,Pr_a] = Air_prop(Ta_avg,Pa)
Ta_avg_K=Ta_avg+273.15;
rho_a=0.4354*(800/(Ta_avg_K))*(Pa/1.013); %[kg/m^3] Uses ideal gas law and density at 800K and atm press.
air=csvread('air_data.csv',1,0); %Air data from Incropera
cp_a=interp1(air(:,1),air(:,2),Ta_avg_K); %[kJ/kg*K]
mu_a=interp1(air(:,1),air(:,3),Ta_avg_K); %[N*s/m^2]
k_a=interp1(air(:,1),air(:,4),Ta_avg_K); %[W/m*K]
Pr_a=cp_a*mu_a/k_a; %dimensionless
