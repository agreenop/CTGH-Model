%This program will create a 2-D model (radial and azimuthal) of a section
%of the CTGH assuming the gas is only flowing in the radial direction.
%This program is used in the 3-D model.
function [T_l_out,T_g_out,P_l_out,T_l,T_g,P_l,P_g,Q,epsilon]=CTGH_2D_calculations(m_l_2_D,m_g_2_D)
load('THEEM_3D_Input.mat');
i=1;
[tubes_vol,N_T,N_L,tubes,D_in,L,H,k_t,rho_t,Cp_t,R_curv,loops,spacers,section,bundles]=CTGH_geom(tube_material,D_out,t,ST,SL,entry,i); %Establishes geometry and material of tubes
Q=zeros(loops*entry+spacers,108); %Establish grid size of system
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
Re_l_matrix=zeros(size(Q));
De_l_matrix=zeros(size(Q));
T_l_out=zeros(entry,1); %Matrix of outlet temperatures for liquid
P_l_out=zeros(entry,1); %Matrix of outlet pressures for liquid
m_l_vol=m_l_2_D/(entry);%Mass flow of liquid through all tubes per volume
m_l_t=m_l_vol/(N_T*N_L); %Mass flow of coolant per tube assuming even distribution
m_g_vol=m_g_2_D/(size(Q,2)); %Mass flow of gas per volume 
for j=1:size(T_g,2)
    T_g(1,j)=T_g_in; %Gas inlet temperature at interior of CTGH
    P_g(1,j)=P_g_in; %Gas inlet pressure at interior of CTGH
end
BC_g=nnz(T_g);%Number of 'cool' gas entry points into CTGH
for j=1:round(size(Q,2)/entry):size(Q,2)
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
for j_entry=1:round(size(Q,2)/entry):size(Q,2) %Equally distributes entry points
 i=size(Q,1); %Starts loop at outer row
 count_move=2; %Counts spaces until next entry point
while i>0
    for j0=j_entry:size(Q,2)+j_entry-1
%This section allows for the spiral motion of the CTGH
        [tubes_vol,N_T,N_L,tubes,D_in,L,H,k_t,rho_t,Cp_t,R_curv,loops,spacers,section,bundles]=CTGH_geom(tube_material,D_out,t,ST,SL,entry,i); 
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
            if i==size(Q,1)-3 || i==size(Q,1)-8
                i1=i-2;
            end
         else
             i1=i;
         end
         if i1==0
            i1=1;
         end
%This section assigns the values of the coefficient matrix that will solve for the temperatures and heat transfer for each volume.        
       l1=(i-1)*size(T_l,2)+j; %Placement of T_l(i,j) coefficient
       l2=(i1-1)*size(T_l,2)+j1; %Placement of T_l(i1,j1) coefficient
       g1=numel(T_l)+(i-1)*size(T_g,2)+j; %Placement of T_g(i,j) coefficient
       g2=numel(T_l)+(i)*size(T_g,2)+j; %Placement of T_g(i+1,j) coefficient
       q1=numel(T_l)+numel(T_g)+(i-1)*size(Q,2)+j; %Placement of Q coefficient
       %EQ1: 0=m_l_vol*Cp_l*(T_l(i,j)-T_l(i,j+1))-Q(i,j);
       count=count+1; %Tracks number of equations
       [UA,Cp_l,Cp_g,mu_l,rho_l,u_max_app,~,Re_g,h_g,Area,Re_l,f_l,De_l]=heat_properties(inlet_prop,gas,liquid,tube_material,D_out,t,ST,SL,T_l_in,T_g_in,P_l_in,P_g_in,T_g,T_l,P_g,P_l,m_g_vol,i,j,i1,j1,m_l_t,model_selection,entry);
       UA_matrix(i,j)=UA; %Records UA values for this volume
       U_matrix(i,j)=UA/Area; %Records U values for this volume
       A_matrix(i,j)=Area; %Records surface area for each volume
       Velocity_matrix(i,j)=u_max_app;
       Re_g_matrix(i,j)=Re_g;
       h_g_matrix(i,j)=h_g;
       Re_l_matrix(i,j)=Re_l;
       De_l_matrix(i,j)=De_l;
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
       if i==size(T_g,1)-5||i==size(T_g,1)-10
          if inlet_prop>1
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
for i=[size(T_g,1)-5,size(T_g,1)-10]
    for j=1:size(T_g,2)
            g1=numel(T_l)+(i-1)*size(T_g,2)+j; %Placement of T_g(i,j) coefficient
            g2=numel(T_l)+(i)*size(T_g,2)+j; %Placement of T_g(i+1,j) coefficient
            count=count+1;
        if inlet_prop==1
            A(count,g1)=1;
            A(count,g2)=-1;
        elseif inlet_prop>1
                if j>=5 && j<=13
                    T_mix=mean(T_g(i,5:13));
                    A(count,g2)=1;
                    B(count)=T_mix;
                elseif j>=14 && j<=22
                    T_mix=mean(T_g(i,14:22));                    
                    A(count,g2)=1;
                    B(count)=T_mix;
                elseif j>=23 && j<=31
                    T_mix=mean(T_g(i,23:31));
                    A(count,g2)=1;
                    B(count)=T_mix;
                elseif j>=32 && j<=40
                    T_mix=mean(T_g(i,32:40));
                    A(count,g2)=1;
                    B(count)=T_mix;
                elseif j>=41 && j<=49
                    T_mix=mean(T_g(i,41:49));
                    A(count,g2)=1;
                    B(count)=T_mix;
                elseif j>=50 && j<=58
                    T_mix=mean(T_g(i,50:58));
                    A(count,g2)=1;
                    B(count)=T_mix;
                elseif j>=59 && j<=67
                    T_mix=mean(T_g(i,59:67));
                    A(count,g2)=1;
                    B(count)=T_mix;
                elseif j>=68 && j<=76
                    T_mix=mean(T_g(i,68:76));
                    A(count,g2)=1;
                    B(count)=T_mix;
                elseif j>=77 && j<=85
                    T_mix=mean(T_g(i,77:85));
                    A(count,g2)=1;
                    B(count)=T_mix; 
                elseif j>=86 && j<=94
                    T_mix=mean(T_g(i,86:94));
                    A(count,g2)=1;
                    B(count)=T_mix;
                elseif j>=95 && j<=103
                    T_mix=mean(T_g(i,77:85));
                    A(count,g2)=1;
                    B(count)=T_mix;
                elseif j>=1 && j<=4 
                    T_mix=mean([T_g(i,104:108),T_g(i,1:4)]);
                    A(count,g2)=1;
                    B(count)=T_mix;
                elseif j>=104 && j<=108
                    T_mix=mean([T_g(i,104:108),T_g(i,1:4)]);
                    A(count,g2)=1;
                    B(count)=T_mix;
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
              if i==size(T_g,1)-5||i==size(T_g,1)-10
                  P_g(i+1,j)=P_g(i,j);
              else
                %Calculates gas pressure drop across bank of tubes.  See Eq. 7.61 in Incopera 5th Ed.
                [~,~,~,mu_l,rho_l,u_max_app,rho_g]=heat_properties(inlet_prop,gas,liquid,tube_material,D_out,t,ST,SL,T_l_in,T_g_in,P_l_in,P_g_in,T_g,T_l,P_g,P_l,m_g_vol,i,j,i1,j1,m_l_t,model_selection,entry);
                chi=1.15;
                f=0.2;
                P_g(i+1,j)=P_g(i,j)-(N_L*chi*rho_g*f*u_max_app^2/2)*10^-5;
              end
          end
         end
     end
%This condition is met when all of the new temperatures of the liquid based 
%on the new temperature-dependent properties are within 0.001 of all of the 
%corresponding old temperatures based on the old temperature-dependent
%properties.
if mean(mean(abs((T_l-T_l_out_old)./T_l)<.01))==1||inlet_prop==5
    break 
end
T_l_out_old=T_l; %Current matrix of liquid temperature is now old matrix of liquid temperatures
inlet_prop=inlet_prop+1; 
end
Q_actual=sum(sum(Q)); %Sum up heat transfer in all volumes to get total heat transfer
Q_avg=mean(mean(Q(Q~=0)));
UA_avg=mean(mean(UA_matrix(UA_matrix~=0)));
Q(Q==0)=NaN; %Places "NaN" in blank spots in Q matrix
T_l(T_l==0)=NaN; %Places "NaN" in blank spots in T_l matrix
UA_matrix(UA_matrix==0)=NaN;
Re_g_matrix(Re_g_matrix==0)=NaN;
Re_l_matrix(Re_l_matrix==0)=NaN;
h_g_matrix(h_g_matrix==0)=NaN;
U_avg=sum(sum(U_matrix))/nnz(U_matrix);
A_total=sum(sum(A_matrix));
U_matrix(U_matrix==0)=NaN;
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
    P_l_out(entry_number,1)=P_l(1,T_l_out(entry_number,1)); %Records outlet temperature of each loop
    T_l_out(entry_number,1)=T_l(1,T_l_out(entry_number,1)); %Records outlet temperature of each loop
end
T_g_out=T_g(size(T_g,1),:);
T_g_avg_out=mean(T_g(size(T_g,1),:)); %Average gas outlet temperature
% T_g_avg_total=(T_g_avg_out+T_g_in)/2; %Average gas temperature across CTGH
T_l_avg_out=mean(T_l_out); %Avergae liquid outlet temperature
% T_l_avg_total=(T_l_avg_out+T_l_in)/2; %Average liquid temperature across CTGH
UA_total=U_avg*A_total;
LMTD=((T_l_in-T_g_avg_out)-(T_l_avg_out-T_g_in))/log((T_l_in-T_g_avg_out)/(T_l_avg_out-T_g_in));
Q_m=UA_total*LMTD;
% inlet_prop=1;
% [UA,Cp_l,Cp_g]=heat_properties(inlet_prop,gas,liquid,tube_material,D_out,t,ST,SL,T_l_avg_total,T_g_avg_total,P_l_in,P_g_in,T_g,T_l,P_g,P_l,m_g,i,j,i1,j1,m_l_t,model_selection);
% C_min=min(m_g_2_D*Cp_g,m_l_2_D*Cp_l);
% Q_max=C_min*(T_l_in-T_g_in);
% e1=Q_actual/Q_max;
epsilon=Q_actual/Q_m;
% fprintf('The effectiveness of this heat exchanger is %4.4f.\n',epsilon)
% CTGH_plot(T_l,T_g,Q,P_l,P_g,UA_matrix,Re_g_matrix,h_g_matrix,Re_l_matrix,U_matrix,gas,liquid) %Plots the values
% save('THEEM_Output.mat');
% load('THEEM_Output.mat');