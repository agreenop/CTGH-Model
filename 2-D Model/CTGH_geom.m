%This function establishes the geomtry of the CTGH bundle and of the control volume.  This function
%will make it easy to redesign the geomerty if necessary.
function [L,R_curv,H,tubes_vol,N_T,N_L,tubes,D_in,section,L_tube_avg,vol_cells_gap,slice_total,slice_holder,R_co,vol_wid]=CTGH_geom(THEEM_model,i)
%% Input CTGH Geometry Parameters
if strcmp(THEEM_model, '3D')
    load('THEEM_Input_3D.mat');
elseif strcmp(THEEM_model, '2D')
    load('THEEM_Input_2D.mat');
else
    load('THEEM_Input_temp_2D.mat');
end
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
N_L=2*tube_layer; %Number of columns in longitudinal/radial direction per volume cell
vol_wid=N_L*SL*D_out; %Width of a volume element in radial direction
tube_vert_list=divisors(layer_num); %All possible number of tube rows in a volume cell
H_possible=((tube_vert_list)/2)*D_out*ST; %Height of each volume cell based on the number of possible tube rows
[~, min_pos1]=min(abs((H_possible-vol_wid)/vol_wid)); %Determines how close each possible height value is to the width of the volume cell. Making these values as close as possible makes the cell more cubic.
N_T=tube_vert_list(min_pos1); %Average number of rows in transverse/vertical direction per volume cell. 
section=layer_num/N_T; %Number of layers of volume cells vertically per sub-bundle
H=(N_T/2)*D_out*ST; %Height of control volume that includes gap at top of each volume
tubes_vol=N_T*(tube_layer-heat_rod); %Number of tubes in finite volume
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