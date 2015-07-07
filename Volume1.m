%%This program is to model the heat transfer for a finite element of the
%%CTAH.  It will focus on the geometry of the most basic geometry.  Look at 
%%0-D model for CTAH in excel.
%%
%%This will have dummy input variables for now.  However, output of one
%%geomtry will eventually equal input for another geometry.
%%
%% Andrew Greenop 8 June 2015.
clc;clear;
Ts(1)=700; %Fluid (Salt) Inlet temp. in celsius
Ta(1)=418.6; %Gas (Air) Inlet temp. in celsius after going through compressor
Pa(1)=18.76; %Gas inlet pressure in bar after compressor
m_a=418.5; %Assumed air flow through system in kg/s
m_s=480.2; %Assumed salt flow through system assuming 700 deg inlet and 600 deg outlet in kg/s
tubes=13680; %Total number of tubes in CTAH
tubes_vol=23; %Number of tubes in finite volume
N_T=5;
N_L=5;
m_st=m_s/tubes; %Mass flow per tube
Ts_avg=Ts(1);
D_out=0.25*0.0254;
t=.035*.0254;
D_in=D_out-2*t;
L=1.585*0.0254;
H=1.26*0.0254;
%m_a_vol=m_a*(1/36)*(1/12)*(1/9)*(1/4); %Air mass flow through volume. 36 sub bundles that was cut into 12 sections. 4 sections of 5 rows each cut into 9 sections each.
m_a_vol=m_a/36;
Nu_s=3.66; %Nusselt number for fully developed laminar flow in a pipe
Ta_avg=Ta(1);
while 1
[mu_s,cp_s,k_s,rho_s,nu_s,Pr_s] = Flibe_prop(Ts_avg);
h_s=k_s/(Nu_s*D_in);
Pdrop_s=128*mu_s*L*m_st/(rho_s*pi*D_in^4);
Re_s=4*m_st/(pi*D_in*mu_s);
R_s=1/(tubes_vol*pi*D_in*L*h_s);
k_m=13.40;
R_m=log(D_out/D_in)/(2*pi*k_m*L);
SL=1.45;
ST=1.256;
[rho_a,cp_a,mu_a,k_a,Pr_a] = Air_prop(Ta_avg,Pa(1));
Re_a=D_out*m_a_vol/(mu_a*H*L);
%The next section is from Khan paper
C2=(0.588+0.004*ST)*(0.858+0.04*ST-0.008*ST^2)^(1/SL);
Nu_df_a=C2*Re_a^(1/2)*Pr_a^(1/3)+0.001*Re_a; 
C1=(1.21+1.64*N_L^1.44)/(1.87+N_L^1.44);
Nu_a=C1*Nu_df_a;
%End of Khan Paper equations
h_a=k_a*Nu_a/D_out;
R_a=1/(tubes_vol*pi*D_out*L*h_a);
UA=1/(R_s+R_m+R_a);
%This next section will use NTU-epislon method.
C_hot=(m_st*tubes_vol)*cp_s;
C_cold=m_a_vol*cp_a;
C_min=min(C_hot,C_cold);
C_max=max(C_hot,C_cold);
C_r=C_min/C_max;
NTU=UA/C_min;
if C_max == C_cold
    epsilon=(1/C_r)*(1-exp(-C_r*(1-exp(-NTU))));
else
    epsilon=1-exp(-C_r^-1*(1-exp(-C_r*NTU)));
end
Q_max=C_min*(Ts(1)-Ta(1));
Q=epsilon*Q_max;
Ts_out=Ts(1)-Q/C_hot;
Ta_out=Ta(1)+Q/C_cold;
F=1;
LMTD=F*((Ts_out-Ta(1))-(Ts(1)-Ta_out))/log((Ts_out-Ta(1))/(Ts(1)-Ta_out));
AMTD=((Ts_out+Ts(1))/2-(Ta_out+Ta(1))/2);
Q1=UA*LMTD;
Ts_avg_new=(Ts(1)+Ts_out)/2;
Ta_avg_new=(Ta(1)+Ta_out)/2;
if abs(Ta_avg_new-Ta_avg)<0.001
    break
end
Ts_avg=Ts_avg_new;
Ta_avg=Ta_avg_new;
end
