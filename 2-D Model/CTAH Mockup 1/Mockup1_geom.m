%This function establishes the geomtry of the CTGH bundle and of the control volume.  This function
%will make it easy to redesign the geomerty if necessary.
function [tubes_vol,N_T,N_L,tubes,D_in,L,H,k_t,rho_t,Cp_t,R_curv,loops]=Mockup1_geom(tube_material,D_out,t,i)
tubes_vol=23; %Number of tubes in finite volume
N_T=10; %Number of rows in transversal direction per volume cell
N_L=2; %Number of rows in longitudinal direction per volume cell
tubes=40; %Total number of tubes in CTGH
D_in=D_out-2*t; %Inside diameter [m]
loops=4; %Number of times mock-up tubes loop
R_ci=17.5*.0254/2; %Inside radius of coiled bundle [m]
R_co=29.0*.0254/2; %Outside radius of coiled bundle [m]
R_curv=R_ci+(R_co-R_ci)*(i-1)/7; %Radius of curvature of volume element [m]
L=2*pi*R_curv/50; %Length of control volume in direction of coolant flow for 50 rows [m]
H=4*0.0254; %Height of control volume [m]
switch tube_material
    case '316 Stainless Steel'
        k_t=13.40; %316 SS thermal conductivity [W/m*K]
        rho_t=8238; %316 SS density [kg/m^3]
        Cp_t=468; %316 SS specific heat [J/kg*K]
end