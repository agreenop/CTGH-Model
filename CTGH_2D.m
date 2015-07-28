%This program will create a 2-D model (radial and azimuthal) of the CTGH
%assuming the gas is only flowing in the radial direction.  For grid layout,
%see Solidworks model.
%Andrew Greenop June 22, 2015
clc;clear;
[gas,liquid,T_g_in,T_l_in,P_g_in,P_l_in,m_g,m_l,tube_material,D_out,t,SL,ST,entry,program]=CTGH_2D_GUI;
if isequal(program,'Cancel')
    return
end
[tubes_vol,N_T,N_L,tubes,D_in,L,H,k_t,rho_t,Cp_t]=CTGH_geom(tube_material,D_out,t); %Establishes geometry and material of tubes
m_l_2_D=m_l/(36*4);
m_l_t=m_l/tubes; %Mass flow of coolant per tube assuming even distribution
m_l_vol=m_l_2_D/(entry);%m_l_t*tubes_vol; %Mass flow of liquid through all tubes per volume %m_l_t*tubes_vol;
m_g_2_D=m_g;%/(36*4);
m_g_vol=m_g_2_D/(108); %Mass flow of gas per volume 0.333722;
Q=zeros(3*entry,108); %Establish grid size of system
T_l=zeros(size(Q,1),size(Q,2));
T_g=zeros(size(Q,1)+1,size(Q,2));
P_g=zeros(size(T_g));
P_l=zeros(size(T_l));
UA_matrix=zeros(size(Q));
T_l_out=zeros(entry,1);
for j=1:size(T_g,2)
    T_g(1,j)=T_g_in; %Gas inlet temperature at interior of CTAH
    P_g(1,j)=P_g_in;
end
BC_g=nnz(T_g);%Number of 'cool' gas entry points into CTAH
for j=1:round(size(Q,2)/entry):size(Q,2)
    T_l(size(T_l,1),j)=T_l_in; %Coolant inlet temperatures at manifold entry points
    P_l(size(T_l,1),j)=P_l_in;
end
BC_l=nnz(T_l); %Number of 'hot' coolant entry points into CTAH
inlet_prop=1;
T_l_out_old=zeros(size(T_l,1),size(T_l,2));
while (1) %Starts process assuming constant heat transfer properties throughout system. Repeats using properties based on average temperature and pressure of volume cell.
A=zeros(numel(T_l)+numel(T_g)+numel(Q)-BC_l-BC_g,numel(T_l)+numel(T_g)+numel(Q));
B=zeros(size(A,1),1);
count=0;
entry_number=1;
for j_entry=1:round(size(Q,2)/entry):size(Q,2)
 i=size(Q,1);
 count_move=3;
while i>0
    for j0=j_entry:size(Q,2)+j_entry-1
%This section allows for the spiral motion of the CTGH
         count_move=count_move+1;
         if j0>size(Q,2)
            j=j0-size(Q,2);
         else
            j=j0;
         end
         if j+1>size(Q,2)
            j1=j+1-size(Q,2);
         else
             j1=j+1;
         end
         if count_move==round(size(Q,2)/entry)||j1==size(T_l,2)-2
         i1=i-1;
         count_move=0;
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
       %EQ1==0=m_l_vol*Cp_l*(T_l(i,j+1)-T_l(i,j))-Q(i,j);
       count=count+1;
       [UA,Cp_l,Cp_g,mu_l,rho_l]=heat_properties(inlet_prop,gas,liquid,tube_material,D_out,t,ST,SL,T_l_in,T_g_in,P_l_in,P_g_in,T_g,T_l,P_g,P_l,m_g_vol,i,j,i1,j1);
       UA_matrix(i,j)=UA;
       A(count,l2)=-m_l_vol*Cp_l;
       A(count,l1)=m_l_vol*Cp_l;
       A(count,q1)=-1;
       if T_l(i1,j1)==T_l_in
           A(count,l2)=0;
           B(count)=T_l(i1,j1)*m_l_vol*Cp_l;
       elseif T_l(i,j)==T_l_in
           A(count,l1)=0;
           B(count)=-T_l(i,j)*m_l_vol*Cp_l;
       end
       %EQ2==0=m_g_vol*Cp_g*(T_g(i,j)-T_g(i+1,j))-Q(i,j);
       count=count+1;
       A(count,g1)=-m_g_vol*Cp_g;
       A(count,g2)=m_g_vol*Cp_g;
       A(count,q1)=-1;
        if T_g(i,j)==T_g_in
           A(count,g1)=0;
           B(count)=T_g(i,j)*m_g_vol*Cp_g;
       elseif T_g(i+1,j)==T_g_in
           A(count,g2)=0;
           B(count)=-T_g(i+1,j)*m_g_vol*Cp_g;
        end
       %EQ3=0=UA*((T_l(i,j+1)+T_l(i,j))/2-(T_g(i+1,j)+T_g(i,j))/2)-Q(i,j);
       count=count+1;
       A(count,l2)=UA/2;
       A(count,l1)=UA/2;
       A(count,g2)=-UA/2;
       A(count,g1)=-UA/2;
       A(count,q1)=-1;
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
       if P_l(i1,j1)~=P_l_in
          P_l(i1,j1)=P_l(i,j)-((128*mu_l*L*m_l_t)/(rho_l*pi*D_in^4))*10^-5; %Pressure drop across volume of coolant assuming laminar flow. 
       end
       if i==1 && j1==size(T_l,2)-2 %Stops loop at T_l(1,size(T_l,2)-2)
           T_l_out(entry_number,1)=j1;
           entry_number=entry_number+1;
           i=0;
            break %Exit for loop. i=0 fullfills while condition.
       elseif j1==j_entry-3 && i==1
           T_l_out(entry_number,1)=j1;
           entry_number=entry_number+1;
           i=0;
           break
       end
       i=i1;
     end
end
end
A_sparse=sparse(A);
B_sparse=sparse(B);
X=mldivide(A,B); %Solves for variables
[T_l,T_g,Q]=assignment(T_l,T_g,Q,X,T_l_in,T_g_in); %Resinserts variables back into their respective matrix
%This next section calculates the pressure drop of the gas
%through each volume element.
i1=1;
     for i=1:size(P_g,1)-1 
         for j=1:size(P_g,2)
          if P_g(i+1,j)~=P_g_in   
            [~,~,~,mu_l,rho_l,u_max_app,rho_g]=heat_properties(inlet_prop,gas,liquid,tube_material,D_out,t,ST,SL,T_l_in,T_g_in,P_l_in,P_g_in,T_g,T_l,P_g,P_l,m_g_vol,i,j,i1,j1);
            chi=1.15;
            f=0.2;
         P_g(i+1,j)=P_g(i,j)-(N_L*chi*rho_g*f*u_max_app^2/2)*10^-5;
          end
         end
     end
if mean(mean(abs(T_l-T_l_out_old)<.001))==1
    break 
end
T_l_out_old=T_l;
inlet_prop=inlet_prop+1;
end
Q_actual=sum(sum(Q));
Q(Q==0)=NaN;
T_l(T_l==0)=NaN;
for i=1:size(T_l,1)
    for j=1:size(T_l,2)
        if isnan(T_l(i,j))==1
            P_l(i,j)=NaN;
            UA_matrix(i,j)=NaN;
        end
    end
end
for i=1:size(T_g,1)
    for j=1:size(T_g,2)
        if T_g(i,j)==0
            T_g(i,j)=T_g(i-1,j);
        end
    end
end
%This next section will calculate the effectiveness of the heat exchanger.
for entry_number=1:size(T_l_out,1)
    T_l_out(entry_number,1)=T_l(1,T_l_out(entry_number,1));
end
T_g_avg_out=mean(T_g(size(T_g,1),:));
T_g_avg_total=(T_g_avg_out+T_g_in)/2;
T_l_avg_out=mean(T_l_out);
T_l_avg_total=(T_l_avg_out+T_l_in)/2;
inlet_prop=1;
[UA,Cp_l,Cp_g]=heat_properties(inlet_prop,gas,liquid,tube_material,D_out,t,ST,SL,T_l_avg_total,T_g_avg_total,P_l_in,P_g_in,T_g,T_l,P_g,P_l,m_g);
C_min=min(m_g_2_D*Cp_g,m_l_2_D*Cp_l);
Q_max=C_min*(T_l_in-T_g_in);
epsilon=Q_actual/Q_max;
fprintf('The effectiveness of this heat exchanger is %4.4f.\n',epsilon)
CTGH_plot(T_l,T_g,Q,P_l,P_g,UA_matrix)