%This function establishes the geomtry of the CTGH bundle and of the control volume.  This function
%will make it easy to redesign the geomerty if necessary.
function [tubes_vol,N_T,N_L,tubes,D_in,L,H,k_t,rho_t,Cp_t,R_curv,loops,spacers,section,bundles]=CTGH_geom(tube_material,D_out,t,ST,SL,entry,i)
%% Volume Cell Parameters
N_T=5; %Number of rows in transversal/vertical direction per volume cell
N_L=5; %Number of columns in longitudinal/radial direction per volume cell
%% General CTGH Parameters
loops=3; %Number of times tube loops around CTGH
row_num=20; %Number of tube rows per sub-bundle
tube_row=5; %Number of tubes per row per manifold pipe
heat_rod=1/2; %Number of heater rods in each tube row (1 every 2 rows)
bundles=36; %Number of sub-bundles in CTGH
spacers=2; %Number of spacer gaps in each sub-bundle (allows for air mixing)
spacer_width=0.038; %Width of each spacer gap based on tie rod diameter [m]
R_ci=1.324/2; %Inside radius of coiled bundle [m]
%% Calculated CTGH Geometry Parameters
tubes_vol=N_T*(tube_row-heat_rod); %Number of tubes in finite volume
vol_wid=SL*D_out+D_out; %Width of a volume element in radial direction
section=row_num/N_T; %Number of rows of volume cells vertically per sub-bundle 
tubes_manifold=row_num*(tube_row-heat_rod); %Number of tubes per manifold per sub-bundle
tubes=entry*tubes_manifold*bundles; %Total number of tubes in CTGH
D_in=D_out-2*t; %Inside diameter [m]
tube_col=entry*(tube_row*2)*loops; %Number of tube columns, inlcuding heating rods & accounting for staggered arrangement
bank_depth=tube_col*SL*D_out+spacers*spacer_width; %Depth of tube bank based on tubes, tube spacing, and spacer gaps [m]
% R_co=R_ci+bank_depth; %Average tube bundle outside radius [m]
R_co=103.99*.0254/2; %Outside radius of coiled bundle [m]
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
% H=2*D_out*ST+D_out; %Height of control volume [m]
H=5*ST*D_out; %Height of control volume [m]
R_curv=R_ci+(R_co-R_ci)*(i-1/2)/13; %Radius of curvature of volume element
% if i<=entry
%     R_curve=R_ci+(vol_wid)*(i-1/2);
% elseif i==entry+1
%     R_ci+(vol_wid)*entry+spacer_width/2;
% elseif i<=2*entry+1
%     R_curve=R_ci+(vol_wid)*entry+spacer_width+(N_L*D_out)*(i-11/2);
% elseif i==2*entry+2
%     R_curve=R_ci+(vol_wid)*entry*2+(3/2)*spacer_width;
% elseif i<=3*entry+2
%     R_curve=R_ci+(vol_wid)*entry*2+2*spacer_width+(N_L*D_out)*(i-21/2);
% end
% L=R_curve*2*pi/108; %Length of control volume in direction of coolant flow [m]
%% Tube Material Properties
switch tube_material
    case '316 Stainless Steel'
        k_t=13.40; %316 SS thermal conductivity [W/m*K]
        rho_t=8238; %316 SS density [kg/m^3]
        Cp_t=468; %316 SS specific heat [J/kg*K]
end