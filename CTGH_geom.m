%This function establishes the geomtry of the CTAH bundle and of the control volume.  This function
%will make it easy to redesign the geomerty if necessary.
function [tubes_vol,N_T,N_L,tubes,D_in,L,H,k_t,rho_t,Cp_t]=CTGH_geom(tube_material,D_out,t)
tubes_vol=23; %Number of tubes in finite volume
N_T=5; %Number of rows in transversal direction
N_L=5; %Number of rows in longitudinal direction
tubes=13680; %Total number of tubes in CTAH
D_in=D_out-2*t; %Inside diameter [m]
L=1.585*0.0254; %Length of control volume in direction of coolant flow [m]
H=1.26*0.0254; %Height of control volume [m]
switch tube_material
    case '316 Stainless Steel'
        k_t=13.40; %316 SS thermal conductivity [W/m*K]
        rho_t=8238; %316 SS density [kg/m^3]
        Cp_t=468; %316 SS specific heat [J/kg*K]
end