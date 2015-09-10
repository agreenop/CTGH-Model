%This function establishes the geomtry of the CTGH bundle and of the control volume.  This function
%will make it easy to redesign the geomerty if necessary.
function [tubes_vol,N_T,N_L,tubes,D_in,L,H,k_t,rho_t,Cp_t]=CTGH_geom(tube_material,D_out,t,i)
tubes_vol=23; %Number of tubes in finite volume
N_T=5; %Number of rows in transversal direction per volume cell
N_L=5; %Number of rows in longitudinal direction per volume cell
tubes=13680; %Total number of tubes in CTGH
D_in=D_out-2*t; %Inside diameter [m]
if i==1    
    L=1.585*0.0254; %Length of control volume in direction of coolant flow [m]
elseif i==2
    L=1.683*0.0254;
elseif i==3
    L=1.788*0.0254;
elseif i==4
    L=1.894*0.0254;
elseif i==5||i==6
    L=2.076*0.0254;
elseif i==7
    L=2.181*0.0254;
elseif i==8
    L=2.286*0.0254;
elseif i==9||i==10
    L=2.392*0.0254;
elseif i==11
    L=2.574*0.0254;
elseif i==12
    L=2.679*0.0254;
elseif i==13
    L=2.784*0.0254;
elseif i>=14
    L=2.89*0.0254;
end
H=1.26*0.0254; %Height of control volume [m]
switch tube_material
    case '316 Stainless Steel'
        k_t=13.40; %316 SS thermal conductivity [W/m*K]
        rho_t=8238; %316 SS density [kg/m^3]
        Cp_t=468; %316 SS specific heat [J/kg*K]
end