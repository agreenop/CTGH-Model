%This function will calculate the thermophysical and maximum allowable
%stress for 316 Stainless Steel.  The is based on ASME BPVC Section II Part
%D. Temperature is in degrees C.
function [k_t,rho_t,Cp_t,S_allow] = SS_316_prop(T_max)
k_t=13.40; %316 SS thermal conductivity [W/m*K]
rho_t=8238; %316 SS density [kg/m^3]
Cp_t=468; %316 SS specific heat [J/kg*K]
if exist('T_max','var')==1
    stress=csvread('316SS_max_stress_data.csv',1,0);
    S_allow=10*interp1(stress(:,1),stress(:,2),T_max); %Max allowable stress in bar
end