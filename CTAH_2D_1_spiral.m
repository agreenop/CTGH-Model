%This program will create a 2-D model (radial and azimuthal) of the CTGH
%assuming the gas is only flowing in the radial direction.  For grid layout,
%see Solidworks model.
%Andrew Greenop June 22, 2015
clc;clear;
[T_g_in,T_l_in,P_g_in,P_l_in,m_g,m_l]=inlet_cond;
[tubes_vol,N_T,N_L,tubes,D_out,D_in,L,H,SL,ST,k_t,rho_t,Cp_t]=CTGH_geom; %Establishes geometry and material of tubes
m_l_t=m_l/tubes; %Mass flow of coolant per tube assuming even distribution
m_l_vol=m_l_t*tubes_vol; %Mass flow of liquid through all tubes per volume
m_g_vol=m_g/(108); %Mass flow of gas per volume
Q=zeros(12,108); %Establish grid size of system
T_l=zeros(size(Q,1),size(Q,2));
T_g=zeros(size(Q,1)+1,size(Q,2));
P_g=zeros(size(T_g));
P_l=zeros(size(T_l));
for j=1:size(T_g,2)
    T_g(1,j)=T_g_in; %Gas inlet temperature at interior of CTAH
    P_g(1,j)=P_g_in;
end
BC_g=nnz(T_g);%Number of 'cool' gas entry points into CTAH
for j=1%,28,55,82]
    T_l(size(T_l,1),j)=T_l_in; %Coolant inlet temperatures at manifold entry points
    P_l(size(T_l,1),j)=P_l_in;
end
BC_l=nnz(T_l); %Number of 'hot' coolant entry points into CTAH
inlet_prop=1;
[UA,Cp_l,Cp_g]=heat_properties(inlet_prop,T_l_in,T_g_in,P_l_in,P_g_in,T_g,T_l,P_g,P_l,m_g_vol);
T_g(2,size(Q,2))=T_g(1,size(Q,2)); %No heat exchanger at grid spacing (1,size(Q,2)).
T_l_out_old=T_l_in;
while (1) %Starts process assuming constant heat transfer properties throughout system. Repeats using properties based on temperature and pressure of volume cell.
A=zeros(numel(T_l)+numel(T_g)+numel(Q)-BC_l-BC_g,numel(T_l)+numel(T_g)+numel(Q));
B=zeros(size(A,1),1);
count=0;
i=size(Q,1);
while i>0
    for j=1:size(Q,2)
%This section allows for the spiral motion of the CTGH       
        if j+1>size(Q,2)
            j1=j+1-size(Q,2);
            i1=i-1;
        else
        j1=j+1;
        i1=i;
        end
        if i1==0 %Stops loop at T_l(1,size(T_l,2))
            i=0; 
            break %Exit for loop. i=0 fullfills while condition.
        end
%This section calculates the temperature and heat transfer for each volume.        
       l1=(i-1)*size(T_l,2)+j; %Placement of T_l(i,j) coefficient
       l2=(i1-1)*size(T_l,2)+j1; %Placement of T_l(i1,j1) coefficient
       g1=numel(T_l)+(i-1)*size(T_g,2)+j; %Placement of T_g(i,j) coefficient
       g2=numel(T_l)+(i)*size(T_g,2)+j; %Placement of T_g(i+1,j) coefficient
       q1=numel(T_l)+numel(T_g)+(i-1)*size(Q,2)+j; %Placement of Q coefficient
       %EQ1==0=m_l_vol*Cp_l*(T_l(i,j+1)-T_l(i,j))-Q(i,j);
       count=count+1;
       [UA,Cp_l,Cp_g,mu_l,rho_l]=heat_properties(inlet_prop,T_l_in,T_g_in,P_l_in,P_g_in,T_g,T_l,P_g,P_l,m_g_vol,i,j,i1,j1);
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
           B(count)=T_g(i,j)*UA/2+B(count);
       elseif T_g(i+1,j)==T_g_in
           A(count,g2)=0;
           B(count)=T_g(i+1,j)*UA/2+B(count);
       end 
       if P_l(i1,j1)~=P_l_in
          P_l(i1,j1)=P_l(i,j)-((128*mu_l*L*m_l_t)/(rho_l*pi*D_in^4))*10^-5; %Pressure drop across volume of coolant assuming laminar flow. 
       end
    end
   i=i-1;
end
A1=A(:,any(A));
X=linsolve(A1,B); %Solves for variables
[T_l,T_g,Q]=assignment(T_l,T_g,Q,X,T_l_in,T_g_in); %Resinserts variables back into their respective matrix
%This next section calculates the pressure drop of the gas
%through each volume element.
i1=1;
     for i=1:size(P_g,1)-1 
         for j=1:size(P_g,2)
          if P_g(i+1,j)~=P_g_in   
            [UA,Cp_l,Cp_g,mu_l,rho_l,u_max_app,rho_g]=heat_properties(inlet_prop,T_l_in,T_g_in,P_l_in,P_g_in,T_g,T_l,P_g,P_l,m_g_vol,i,j,i1,j1);
            chi=1.15;
            f=0.2;
         P_g(i+1,j)=P_g(i,j)-(N_L*chi*rho_g*f*u_max_app^2/2)*10^-5;
          end
         end
     end
if abs(T_l(1,size(T_l,2))-T_l_out_old)<0.000001
   break 
end
T_l_out_old=T_l(1,size(T_l,2));
inlet_prop=inlet_prop+1;
end
% disp(T_l);
% disp(T_g);
% disp(Q);
% disp(P_l);
% disp(P_g);
CTGH_plot(T_l,T_g,Q,P_l,P_g)