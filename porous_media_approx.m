%This function determines the mass flow rate distribution of the gas
%flowing through each bundle of the CTGH bundle.  It is based off the idea
%of a porous media approximation.  This is assuming the permeability of the
%azimuthal flow is a fraction of the radial flow, which is the primary
%direction of flow (v=k*u where k<1).  For simplicity the flow is assumed
%isothermal for distrubtion.
function m_g_vol=porous_media_approx(m_g_2_D,T_g,slice_holder,slice_total,gaps_position,THEEM_model)
if strcmp(THEEM_model, '3D')
    load('THEEM_Input_3D.mat');
elseif strcmp(THEEM_model, '2D')
    load('THEEM_Input_2D.mat');
else
    load('THEEM_Input_temp_2D.mat');
end
k_perm=0.01; %Ratio of azimuthal permeability to radial permeability
% k_perm_gap=10; %Ratio of permeability of gap to tube bundle
m_g_vol=zeros(size(T_g)); %Initialize mass flow matrix
m_g_vol(1,:)=m_g_2_D/(size(m_g_vol,2)); %Evenly distribute flow initially between each volume at the inside of the annulus
if mod(slice_holder,2)==1 %If odd number of volume cells between each tube holder
    gap_initial=(slice_holder+1)/2;
else %If even number of volume cells
    gap_initial=slice_holder/2+1;
end
gap_count=1;
gap_matrix=zeros(slice_holder,tube_holders); %Initialize 
for gap_index=gap_initial:slice_total+gap_initial-1
  	if gap_index>slice_total
       	gap_pos_az=gap_index-slice_total;
    else
        gap_pos_az=gap_index;
    end
    gap_matrix(gap_count)=gap_pos_az;
    gap_count=gap_count+1;
end
slice_initial_pos=gap_initial; %The initial volume, which will be the middle of each slice, is based on the offset of the tube holders in relation to the first manifold.
for slice=1:tube_holders %Look at each slice.  (Tube holders prevent gas from mixing between slices)
    for i=2:size(m_g_vol,1)
        if any(i==gaps_position)
            vol_collect=zeros(size(gap_matrix,1),1);
            for gap_collect=1:size(gap_matrix,1)
                vol_collect(gap_collect)=m_g_vol(gap_matrix(gap_collect,slice),i);
            end
            m_g_avg=mean(vol_collect);
            for j=gap_matrix(1,slice):gap_matrix(size(gap_matrix,1),slice)
            m_g_vol(i,j)=m_g_avg;
            end
        end
        for slice_pos=slice_initial_pos:-1:1 %Start at middle and work backwards in each slice
            j=gap_matrix(slice_pos,slice); %Azimuthal position determined by gap matrix
            j1=gap_matrix(slice_pos+1,slice);
            if slice_pos==slice_initial_pos %Need to split initial (middle) volume in half
                m_g_vol(i,j)=(1-k_perm)*m_g_vol(i-1,j);
            elseif slice_pos==slice_initial_pos-1
                m_g_vol(i,j)=(1-k_perm)*(m_g_vol(i-1,j)+k_perm*0.5*m_g_vol(i,j1));
            elseif slice_pos==1
                m_g_vol(i,j)=(m_g_vol(i-1,j)+k_perm*m_g_vol(i,j1));
            else
                m_g_vol(i,j)=(1-k_perm)*(m_g_vol(i-1,j)+k_perm*m_g_vol(i,j1));
            end
        end
        for slice_pos=slice_initial_pos+1:slice_holder %Start at middle and work forwards in each slice
            j=gap_matrix(slice_pos,slice); %Azimuthal position determined by gap matrix
            j1=gap_matrix(slice_pos-1,slice);
            if slice_pos==slice_initial_pos %Need to split initial (middle) volume in half
                m_g_vol(i,j)=(1-k_perm)*m_g_vol(i-1,j);
            elseif slice_pos==slice_initial_pos+1
                m_g_vol(i,j)=(1-k_perm)*(m_g_vol(i-1,j)+k_perm*0.5*m_g_vol(i,j1));
            elseif slice_pos==slice_holder
                m_g_vol(i,j)=(m_g_vol(i-1,j)+k_perm*m_g_vol(i,j1));
            else
                m_g_vol(i,j)=(1-k_perm)*(m_g_vol(i-1,j)+k_perm*m_g_vol(i,j1));
            end
        end
    end
end