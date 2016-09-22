%This function establishes the geomtry of the CTGH bundle and of the control volume.  This function
%will make it easy to redesign the geomerty if necessary.
function [L,R_curv,H,tubes_vol,N_T,N_L,tubes,D_in,k_t,rho_t,Cp_t,section,L_tube_avg,vol_cells_gap,slice_total,slice_holder,R_co,vol_wid]=CTGH_geom(THEEM_model,i)
%% General CTGH Parameters
if strcmp(THEEM_model, '3D')
    load('THEEM_Input_3D.mat');
else
    load('THEEM_Input_2D.mat');
end
%loops=3; %Number of times tube loops around CTGH
%layer_num=40; %Number of tube layers per sub-bundle
%tube_layer=3; %Number of tubes per layer per manifold pipe
heat_rod=1/2; %Number of heater rods in each tube layer (1 every 2 layer)
% bundles=36; %Number of sub-bundles in CTGH
% spacers=2; %Number of spacer gaps from tie rods in each sub-bundle (allows for air mixing)
% spacer_width=0.038; %Width of each spacer gap based on tie rod diameter [m]
% tube_holders=12; %Number of tube holders per sub-bundle
% R_ci=1.324/2; %Inside radius of coiled bundle [m]
%% Calculated CTGH Geometry Parameters
tubes_manifold=layer_num*(tube_layer-heat_rod); %Number of tubes per manifold per sub-bundle
tubes=entry*tubes_manifold*bundles; %Total number of tubes in CTGH
D_in=D_out-2*t; %Inside diameter [m]
tube_col=entry*(tube_layer*2)*loops; %Number of tube columns, inlcuding heating rods & accounting for staggered arrangement
bank_depth=tube_col*SL*D_out+spacers*spacer_width; %Depth of tube bank based on tubes, tube spacing, and spacer gaps [m]
R_co=R_ci+bank_depth; %Average tube bundle outside radius [m]
R_c_avg=R_ci+0.5*bank_depth; %Radius to middle of tube bundle [m]
L_tube_avg=loops*pi*2*R_c_avg; %Average length of each tube in bundle [m]
vol_cells_gap=ceil(loops*entry/(spacers+1)); %Number of volume cells between each tie rod gap
%% General Volume Cell Parameters
N_L=tube_layer; %Number of columns in longitudinal/radial direction per volume cell
vol_wid=2*(N_L-1/2)*SL*D_out+D_out; %Width of a volume element in radial direction
tube_vert_list=divisors(layer_num); %All possible number of tube rows in a volume cell
H_possible=((tube_vert_list-1)/2)*D_out*ST+D_out; %Height of each volume cell based on the number of possible tube rows
[~, min_pos1]=min(abs((H_possible-vol_wid)/vol_wid)); %Determines how close each possible height value is to the width of the volume cell. Making these values as close as possible makes the cell more cubic.
N_T=tube_vert_list(min_pos1); %Number of rows in transverse/vertical direction per volume cell. 
section=layer_num/N_T; %Number of layers of volume cells vertically per sub-bundle
H=((N_T-1)/2)*D_out*ST+D_out; %Height of control volume 
tubes_vol=N_T*(N_L-heat_rod); %Number of tubes in finite volume
slice_holder_list=1:20; %Calculates how many slices are needed between each tube holder
L_possible=(2*pi*R_c_avg)./(tube_holders*slice_holder_list); 
[~, min_pos2]=min(abs((L_possible-vol_wid)/vol_wid)); %Determines how close each possible length value is to the width of the volume cell. Making these values as close as possible makes the cell more cubic.
slice_holder=slice_holder_list(min_pos2); %Number of volume elements between each holder azimuthally
slice_total=slice_holder*tube_holders; %Number of volume elements azimuthally around total bundle
%% Specific Volume Cell Parameters Based on Position
if i<=vol_cells_gap %Any volume element before the 1st tie rod gap
    R_curv=R_ci+(vol_wid)*(i-1/2);
else
    for tr=1:spacers %Include all tie rod gaps
        if i==tr*(vol_cells_gap+1) %Volument element in a tie rod gap
            R_curv=R_ci+(vol_wid)*tr*vol_cells_gap+spacer_width*(tr-1)+spacer_width/2;
            break
        elseif i>tr*(vol_cells_gap+1) && i<(tr+1)*(vol_cells_gap+1) %Volument element after the first tie rod gap
            R_curv=R_ci+vol_wid*tr*vol_cells_gap+spacer_width*tr+vol_wid*(i-(tr*(vol_cells_gap+1)+1/2));
            break
        end
    end
end    
L=R_curv*2*pi/slice_total; %Length of control volume in direction of coolant flow [m]
%% Tube Material Properties
switch tube_material
    case '316 Stainless Steel'
        k_t=13.40; %316 SS thermal conductivity [W/m*K]
        rho_t=8238; %316 SS density [kg/m^3]
        Cp_t=468; %316 SS specific heat [J/kg*K]
end