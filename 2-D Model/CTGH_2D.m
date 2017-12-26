%This program will create a 2-D model (radial and azimuthal) of the CTGH
%assuming the gas is only flowing in the radial direction.  For grid layout,
%see Solidworks model.
%Andrew Greenop June 22, 2015
function [T_l_out,T_g_out,P_l_out,T_l,T_g,P_l,P_g,Q,F_factor,epsilon,U_avg,A_total]=CTGH_2D(THEEM_model,m_l_2_D_exact,m_g_2_D_exact)
if strcmp(THEEM_model, '3D')
    load('THEEM_Input_3D.mat');
elseif strcmp(THEEM_model, '2D')
    load('THEEM_Input_2D.mat');
else
    load('THEEM_Input_temp_2D.mat');
end
i=1;
[L,R_curv,H,tubes_vol,N_T,N_L,tubes,D_in,section,L_tube_avg,vol_cells_gap,slice_total,slice_holder,R_co,vol_wid]=CTGH_geom(THEEM_model,i);%Establishes geometry of tubes
Q=zeros(loops*entry+spacers,slice_total); %Establish grid size of system
T_l=zeros(size(Q,1),size(Q,2));
T_g=zeros(size(Q,1)+1,size(Q,2));
P_g=zeros(size(T_g));
P_l=zeros(size(T_l));
UA_matrix=zeros(size(Q));
U_matrix=zeros(size(Q));
A_matrix=zeros(size(Q));
Velocity_matrix=zeros(size(Q));
Re_g_matrix=zeros(size(Q));
h_g_matrix=zeros(size(Q));
h_l_matrix=zeros(size(Q));
Re_l_matrix=zeros(size(Q));
De_l_matrix=zeros(size(Q));
Nu_l_matrix=zeros(size(Q));
Nu_g_matrix=zeros(size(Q));
i1_matrix=zeros(size(Q));
j1_matrix=zeros(size(Q));
T_s_in_matrix=zeros(size(Q));
T_s_out_matrix=zeros(size(Q));
T_l_out=zeros(entry,1); %Matrix of outlet temperatures for liquid
P_l_out=zeros(entry,1); %Matrix of outlet pressures for liquid
if strcmp(THEEM_model, '3D') %For 3-D calculations, calculates for each bundle
    m_l_2_D=m_l_2_D_exact; %Mass flow split between 36 bundles, which are split into 4 each
    m_g_2_D=m_g_2_D_exact;
else % For 2-D calculations, averages flow rates across bundles
    m_l_2_D=m_l/(bundles*section); %Mass flow split between 36 bundles, which are split into 4 each
    m_g_2_D=m_g/(bundles*section);
end
m_l_vol=m_l_2_D/(entry);%Mass flow of liquid through all tubes per volume
m_l_t=m_l_vol/(tubes_vol); %Mass flow of coolant per tube assuming even distribution
m_g_vol=m_g_2_D/(size(Q,2)); %Mass flow of gas per volume
gaps_position=zeros(1,spacers);
for gap=1:spacers %Start the count for tie rod gaps
    gaps_position(gap)=size(Q,1)-(vol_cells_gap-1)-(gap-1)*(1+vol_cells_gap);%Radial locations of volumes adjacent to tie rod gaps
end
% m_g_vol=porous_media_approx(m_g_2_D,T_g,slice_holder,slice_total,gaps_position,tube_holders); %External function determines isothermal gas flow using a porous media approximation
for j=1:size(T_g,2)
    T_g(1,j)=T_g_in; %Gas inlet temperature at interior of CTGH
    P_g(1,j)=P_g_in; %Gas inlet pressure at interior of CTGH
end
BC_g=nnz(T_g);%Number of 'cool' gas entry points into CTGH
entry_step=round(size(Q,2)/entry);%Establishes the number of azimuthal elements that should be between each liquid manifold. Round forces it to be an integer.
max_entry=size(Q,2)+2-entry_step; %Prevents the code from establishing n+1 entry points.  Add 1 to make the 1st point equal to total+1.  Add another 1 to account for round functioning rounding up instead of down.
for j=1:entry_step:max_entry
    T_l(size(T_l,1),j)=T_l_in; %Coolant inlet temperatures at manifold entry points
    P_l(size(T_l,1),j)=P_l_in; %Coolant inlet pressure at manifold entry points
end
BC_l=nnz(T_l); %Number of 'hot' coolant entry points into CTGH
inlet_prop=1; %Determines which temperature heat properties are taken at for liquid and gas.
T_l_out_old=zeros(size(T_l,1),size(T_l,2));
while (1) %Starts process assuming constant heat transfer properties throughout system. Repeats using properties based on average temperature and pressure of volume cell.
A=zeros(numel(T_l)+numel(T_g)+numel(Q)-BC_l-BC_g,numel(T_l)+numel(T_g)+numel(Q));
B=zeros(size(A,1),1);
count=0; %Count that tracks the number of equations used
entry_number=1; %Tracks which entry point/loop that i currently being calculated
for j_entry=1:entry_step:max_entry %Equally distributes entry points
 i=size(Q,1); %Starts loop at outer row
 count_move=2; %Counts spaces until next entry point
while i>0
    for j0=j_entry:size(Q,2)+j_entry-1
%This section allows for the spiral motion of the CTGH
        [L]=CTGH_geom(THEEM_model,i); 
         count_move=count_move+1;
         if j0>size(Q,2)
            j=j0-size(Q,2); %Resets j to 1 after full rotation around CTGH
         else
            j=j0;
         end
         if j+1>size(Q,2)
            j1=j+1-size(Q,2); %Resets j1 to 1 after full rotation around CTGH
         else
             j1=j+1; %j1 is the next grid space azimuthally from j
         end
         if count_move==round(size(Q,2)/entry)||j==size(T_l,2)-2 
            i1=i-1; %Moves down a row when in line with an entry point
            count_move=0; %Resets spacing from entry point
            if any(i==gaps_position)==1 %If at a volumes adjacent to tie rod gap
                i1=i-2; %Skips over the tie rod gap
            end            
         else
             i1=i;
         end
         if i1==0
            i1=1;
         end
       %These matrices store the position of the next volume for each volume based on the flow path of the liquid. 
       i1_matrix(i,j)=i1; 
       j1_matrix(i,j)=j1;
%This section assigns the values of the coefficient matrix that will solve for the temperatures and heat transfer for each volume.        
       l1=(i-1)*size(T_l,2)+j; %Placement of T_l(i,j) coefficient
       l2=(i1-1)*size(T_l,2)+j1; %Placement of T_l(i1,j1) coefficient
       g1=numel(T_l)+(i-1)*size(T_g,2)+j; %Placement of T_g(i,j) coefficient
       g2=numel(T_l)+(i)*size(T_g,2)+j; %Placement of T_g(i+1,j) coefficient
       q1=numel(T_l)+numel(T_g)+(i-1)*size(Q,2)+j; %Placement of Q coefficient
       %EQ1: 0=m_l_vol*Cp_l*(T_l(i,j)-T_l(i,j+1))-Q(i,j);
       count=count+1; %Tracks number of equations
       [UA,Cp_l,Cp_g,rho_l,u_max_app,rho_g,Re_g,h_g,Area,Re_l,f_l,De_l,h_l,Nu_l,Nu_g]=heat_properties(inlet_prop,T_g,T_l,P_g,m_g_vol,i,j,i1,j1,m_l_t,T_s_out_matrix,T_s_in_matrix,THEEM_model);
       UA_matrix(i,j)=UA; %Records UA values for this volume
       U_matrix(i,j)=UA/Area; %Records U values for this volume
       A_matrix(i,j)=Area; %Records surface area for each volume
       Velocity_matrix(i,j)=u_max_app;
       Re_g_matrix(i,j)=Re_g;
       h_g_matrix(i,j)=h_g;
       h_l_matrix(i,j)=h_l;
       Re_l_matrix(i,j)=Re_l;
       De_l_matrix(i,j)=De_l;
       Nu_g_matrix(i,j)=Nu_g;
       Nu_l_matrix(i,j)=Nu_l;
       A(count,l2)=-m_l_vol*Cp_l;
       A(count,l1)=m_l_vol*Cp_l;
       A(count,q1)=-1;
       %If one of the temperatures is a known value, this section moves it
       %to the solution matrix, B.
       if T_l(i1,j1)==T_l_in
           A(count,l2)=0;
           B(count)=T_l(i1,j1)*m_l_vol*Cp_l;
       elseif T_l(i,j)==T_l_in
           A(count,l1)=0;
           B(count)=-T_l(i,j)*m_l_vol*Cp_l;
       end
       %EQ2: 0=m_g_vol*Cp_g*(T_g(i+1,j)-T_g(i,j))-Q(i,j);
       count=count+1; %Tracks number of equations
       A(count,g1)=-m_g_vol*Cp_g; 
       A(count,g2)=m_g_vol*Cp_g;
       A(count,q1)=-1;
       %If one of the temperatures is a known value, this section moves it
       %to the solution matrix, B.
        if T_g(i,j)==T_g_in
           A(count,g1)=0;
           B(count)=T_g(i,j)*m_g_vol*Cp_g;
       elseif T_g(i+1,j)==T_g_in
           A(count,g2)=0;
           B(count)=-T_g(i+1,j)*m_g_vol*Cp_g;
        end
       %EQ3: 0=UA*((T_l(i,j+1)+T_l(i,j))/2-(T_g(i+1,j)+T_g(i,j))/2)-Q(i,j);
       count=count+1; %Tracks number of equations
       A(count,l2)=UA/2;
       A(count,l1)=UA/2;
       A(count,g2)=-UA/2;
       A(count,g1)=-UA/2;
       A(count,q1)=-1;
       %If one of the temperatures is a known value, this section moves it
       %to the solution matrix, B.
       if T_l(i1,j1)==T_l_in
           A(count,l2)=0;
           B(count)=-T_l(i1,j1)*UA/2;
       elseif T_l(i,j)==T_l_in
           A(count,l1)=0;
           B(count)=-T_l(i,j)*UA/2;
       end
       if T_g(i,j)==T_g_in
           A(count,g1)=0;
           B(count)=T_g(i,j)*UA/2;
       elseif T_g(i+1,j)==T_g_in
           A(count,g2)=0;
           B(count)=T_g(i+1,j)*UA/2;
       end
       if inlet_prop>1
          if any(i+1==gaps_position)==1
            A(count,g1)=0;
            B(count)=T_mix*UA/2;
          end
       end
       %This section calculates the pressure of the liquid in each volume
       %using friction head loss formula.
       if P_l(i1,j1)~=P_l_in
          P_l(i1,j1)=P_l(i,j)-f_l*(L/D_in^5)*(8*m_l_t^2/(pi^2*rho_l))*10^-5; %Pressure drop across volume of coolant. 
       end
       if i==1 && j1==size(T_l,2)-2 %Stops loop at T_l(1,size(T_l,2)-2), a liquid exit point. Only applies to loop with entry point at j=1.
           count=count+1;
           g1=numel(T_l)+(i1-1)*size(T_g,2)+j1; %Placement of final T_g(i1,j1) coefficient in loop
           g2=numel(T_l)+(i1)*size(T_g,2)+j1; %Placement of final T_g(i1+1,j1) coefficient in loop
           q1=numel(T_l)+numel(T_g)+(i1-1)*size(Q,2)+j1; %Placement of final Q(i1,j1) coefficient in loop
           B(count)=m_g_vol*Cp_g*T_g(i1,j1); 
           A(count,g2)=m_g_vol*Cp_g;
           A(count,q1)=-1;
           T_l_out(entry_number,1)=j1; %Records outlet liquid temperature position of this loop in azimuthal direction
           entry_number=entry_number+1; %Moves counter to next loop.
           i=0;
           break %Exit for loop. i=0 fullfills while condition.
       elseif j1==j_entry-3 && i==1 %Stops loop at all other liquid exit points besides loop with entry point at j=1.
           count=count+1;
           g1=numel(T_l)+(i1-1)*size(T_g,2)+j1; %Placement of final T_g(i1,j1) coefficient in loop
           g2=numel(T_l)+(i1)*size(T_g,2)+j1; %Placement of final T_g(i1+1,j1) coefficient in loop
           q1=numel(T_l)+numel(T_g)+(i1-1)*size(Q,2)+j1; %Placement of final Q(i1,j1) coefficient in loop
           B(count)=m_g_vol*Cp_g*T_g(i1,j1);
           A(count,g2)=m_g_vol*Cp_g;
           A(count,q1)=-1;
           T_l_out(entry_number,1)=j1; %Records outlet liquid temperature of this loop
           entry_number=entry_number+1; %Moves counter to next loop.
           i=0;
           break %Exit for loop. i=0 fullfills while condition.
       end
       i=i1;
     end
end
end
for i=gaps_position-1
    for j=1:size(T_g,2)
            g1=numel(T_l)+(i-1)*size(T_g,2)+j; %Placement of T_g(i,j) coefficient
            g2=numel(T_l)+(i)*size(T_g,2)+j; %Placement of T_g(i+1,j) coefficient
            count=count+1;
        if inlet_prop==1
            A(count,g1)=1;
            A(count,g2)=-1;
        elseif inlet_prop>1
            %This section finds how large the tie rod gaps are azimuthally
            %and gives their locations azimuthally.
            if mod(slice_holder,2)==1 %If odd number of volume cells between each tube holder
                gap_initial=(slice_holder+1)/2;
            else %If even number of volume cells
                gap_initial=slice_holder/2+1;
            end
            gap_count=1;
            gap_matrix=zeros(slice_holder,slice_total/slice_holder);
            for gap_index=gap_initial:slice_total+gap_initial-1
            	if gap_index>slice_total
                	gap_pos=gap_index-slice_total;
                else
                	gap_pos=gap_index;
            	end
                        gap_matrix(gap_count)=gap_pos;
                        gap_count=gap_count+1;
            end
            %This loop finds the average gas temperature of each gap 
            for j_gap=1:size(gap_matrix,2)
                if sum(gap_matrix(:,j_gap)==j)~=0 %If j is within the given range
                    T_mix=mean(T_g(i,gap_matrix(:,j_gap)));
                    A(count,g2)=1;
                    B(count)=T_mix;  
                    break
                end
            end
         end
    end
end
X=mldivide(A,B); %Solves for variables (AX=B)
[T_l,T_g,Q]=assignment(T_l,T_g,Q,X,T_l_in,T_g_in); %Resinserts values back into their respective matrix
%This next section calculates the pressure drop of the gas
%through each volume element.
i1=1;
     for i=1:size(P_g,1)-1 
         for j=1:size(P_g,2)
          if P_g(i+1,j)~=P_g_in
              if i==gaps_position-1
                  P_g(i+1,j)=P_g(i,j);
              else
                %Calculates gas pressure drop across bank of tubes.  See Eq. 7.61 in Incopera 5th Ed.
                [~,~,~,~,u_max_app,rho_g,Re_g]=heat_properties(inlet_prop,T_g,T_l,P_g,m_g_vol,i,j,i1,j1,m_l_t,T_s_out_matrix,T_s_in_matrix,THEEM_model);
                [dP_total] = StaggeredPressureDrop(ST,SL,u_max_app,rho_g,N_L,Re_g);
                P_g(i+1,j)=P_g(i,j)-dP_total;
              end
          end
         end
     end
%This section calculates the inside and outside surface temperatures of the
%tubes across the tube bundle.
for i=1:size(Q,1)
    for j=1:size(Q,2)
        if Q(i,j)~=0
           T_s_in_matrix(i,j)=(T_l(i,j)+T_l(i1_matrix(i,j),j1_matrix(i,j)))/2-Q(i,j)/(h_l_matrix(i,j)*A_matrix(i,j)*D_in/D_out);
           T_s_out_matrix(i,j)=(T_g(i,j)+T_g(i+1,j))/2+Q(i,j)/(h_g_matrix(i,j)*A_matrix(i,j));
        end
    end
end
%This condition is met when all of the new temperatures of the liquid based 
%on the new temperature-dependent properties are within 0.1% of all of the 
%corresponding old temperatures based on the old temperature-dependent
%properties.
test_conv=abs((T_l-T_l_out_old)./T_l); %Matrix of percent difference between old temperature values and new values to test for convergence
test_conv(isnan(test_conv))=0; %Eliminates NaN values from dividing by zero
if mean(mean(test_conv<.001))==1||inlet_prop==5
    break 
end
T_l_out_old=T_l; %Current matrix of liquid temperature is now old matrix of liquid temperatures
inlet_prop=inlet_prop+1; 
end
Q_actual=sum(sum(Q)); %Sum up heat transfer in all volumes to get total heat transfer
Q(Q==0)=NaN; %Places "NaN" in blank spots in Q matrix
T_l(T_l==0)=NaN; %Places "NaN" in blank spots in T_l matrix
UA_matrix(UA_matrix==0)=NaN;
Re_g_matrix(Re_g_matrix==0)=NaN;
Re_l_matrix(Re_l_matrix==0)=NaN;
Nu_l_matrix(Nu_l_matrix==0)=NaN;
Nu_g_matrix(Nu_g_matrix==0)=NaN;
h_l_matrix(h_l_matrix==0)=NaN;
h_g_matrix(h_g_matrix==0)=NaN;
U_avg=sum(sum(U_matrix))/nnz(U_matrix);
A_total=sum(sum(A_matrix));
U_matrix(U_matrix==0)=NaN;
T_s_in_matrix(T_s_in_matrix==0)=NaN;
T_s_out_matrix(T_s_out_matrix==0)=NaN;
for i=1:size(T_l,1)
    for j=1:size(T_l,2)
        if isnan(T_l(i,j))==1
            P_l(i,j)=NaN; %Assigns "NaN" to corresponding spots in P_l matrix
            UA_matrix(i,j)=NaN; %Assigns "NaN" to corresponding spots in UA matrix
        end
    end
end
for i=1:size(T_g,1)
    for j=1:size(T_g,2)
        if T_g(i,j)==0
            T_g(i,j)=T_g(i-1,j); %If no heat transfer occurs in a volume, the gas temperature of the previous volume is equal to the current gas temperature.
        end
    end
end
%This next section will calculate the effectiveness of the heat exchanger.
for entry_number=1:size(T_l_out,1)
    P_l_out(entry_number,1)=P_l(1,T_l_out(entry_number,1)); %Records outlet pressure of each loop
    T_l_out(entry_number,1)=T_l(1,T_l_out(entry_number,1)); %Records outlet temperature of each loop
end
%Bundle outlet conditions
T_g_out=T_g(size(T_g,1),:); %Records outlet gas temperatures around the bundle
P_g_out=P_g(size(P_g,1),:); %Records outlet gas pressures around the bundle
T_g_outlet=mean(T_g_out); %Average gas outlet temperature
P_g_outlet=mean(P_g_out); %Average gas outlet pressure
T_g_avg_total=(T_g_outlet+T_g_in)/2; %Average gas temperature across CTGH
P_g_avg_total=(P_g_outlet+P_g_in)/2; %Average gas pressure across CTGH
T_l_outlet=mean(T_l_out); %Average liquid outlet temperature
P_l_outlet=mean(P_l_out); %Average liquid outlet pressures
T_l_avg_total=(T_l_outlet+T_l_in)/2; %Average liquid temperature across CTGH
deltaP_l=P_l_in-P_l_outlet; %Average liquid pressure drop
deltaP_g=P_g_in-P_g_outlet; %Average gas pressure drop
%Bundle effectiveness calculations
UA_total=U_avg*A_total;
LMTD=((T_l_in-T_g_outlet)-(T_l_outlet-T_g_in))/log((T_l_in-T_g_outlet)/(T_l_outlet-T_g_in));
Q_m=UA_total*LMTD; %Heat transfer if perfect counterflow heat exchanger
[Cp_l,Cp_g]=Material_prop(liquid,gas,tube_material,T_l_avg_total,T_g_avg_total,P_g_avg_total);
C_min=min(m_g_2_D*Cp_g,m_l_2_D*Cp_l);
Q_max=C_min*(T_l_in-T_g_in); % Maximum heat transfer possible given inlets 
epsilon=Q_actual/Q_max; %Effectiveness calculation
F_factor=Q_actual/Q_m; %F factor calculation (Comparison with counterflow heat exchanger)
if strcmp(THEEM_model, '2D')
    Q_tot=section*bundles*Q_actual; %Assuming all 2-D cross sections are equal in the bundle, this calculates the overall heat transfer in the bundle.
    fprintf('The effectiveness of this heat exchanger Is %4.4f.\n',epsilon)
    fprintf('The %s outlet temperature is %1.1f%cC.\n',liquid,T_l_outlet,char(176))
    fprintf('The %s outlet temperature is %1.1f%cC.\n',gas,T_g_outlet,char(176))
    fprintf('The %s pressure drop is %1.2f bar.\n',liquid,deltaP_l)
    fprintf('The %s pressure drop is %1.4f bar.\n',gas,deltaP_g)
    fprintf('The estimated overall heat transfer is %1.3e W.\n',Q_tot)
    CTGH_plot(T_l,T_g,Q,P_l,P_g,UA_matrix,Re_g_matrix,h_g_matrix,h_l_matrix,Re_l_matrix,U_matrix,T_s_in_matrix,T_s_out_matrix,gas,liquid,R_ci,vol_cells_gap,vol_wid,spacer_width) %Plots the values
    save('2-D Model/THEEM_Output_2D.mat');
elseif strcmp(THEEM_model, 'Parametric Study')
    Q_tot=section*bundles*Q_actual; %Assuming all 2-D cross sections are equal in the bundle, this calculates the overall heat transfer in the bundle.
    save('Optimization Program/Parametric Study/THEEM_Output_temp_2D.mat')
end