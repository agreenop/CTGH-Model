%Original Author: Jonah Haefner
%Last Modified: 10/27/2015
%Most Reecent Author: Jonah Haefner
%References: Julien Clayton, Lane Carasik, ...
%Modified for THEEM code: Andrew Greenop (12/13/2016)

%%%Units%%%
%ST = dimensionless transverse pitch
%SL = dimensionless longitudinal pitch
%u_max_app = the maximum fluid velocity between tubes in m/s 
%rho_g = density in kg/m^3
%N_L = number of rows in the tube bundle 
%Will work for values of SL = 1.25, 1.5, or 2
%Reynolds number needs to be below 150000 and greater than E2


%coefficients for Euler number (calculated using the power series)
%coefficients come from here: http://www.thermopedia.com/content/1211/#TUBE_BANKS_CROSSFLOW_OVER_FIG2

% These are valid for Reynolds numbers between 2E3 and 2E6 %
% THEEM modifications: Eliminated unneccessary variables, changed to THEEM variables, and added interpolation for exact ratio. 
function [dP_total] = StaggeredPressureDrop(ST,SL,u_max_app,rho_g,N_L,Re_g)
x = double((ST)./(SL));      %Used for correction factor k_1 which I don't have data for.

% Coefficients, c_i, to generate pressure drop coefficients for equilateral triangle banks.
if Re_g < 1E3
c_0 = [.795,.683,.343];      %b = 1.25, 1.5, and 2
c_1 = [.247E3,.111E3,.303E3]; %b = 1.25, 1.5, and 2
c_2 = [.335E3,-.973E2,-.717E5]; %b = 1.25, 1.5, and 2
c_3 = [-.155E4,-.426E3,.88E7]; %b = 1.25, 1.5, and 2
c_4 = [.241E4,.574E3,-.38E9];     %b = 1.25, 1.5, and 2
elseif Re_g >= 1E3 && Re_g < 1E4
c_0 = [.245,.203,.343];      %b = 1.25, 1.5, and 2
c_1 = [.339E4, .248E4, .303E3 ]; %b = 1.25, 1.5, and 2
c_2 = [-.984E7, -.758E7, -.717E5]; %b = 1.25, 1.5, and 2
c_3 = [.132E11, .104E11, .88E7]; %b = 1.25, 1.5, and 2
c_4 = [-.599E13, -.482E13, -.38E9];     %b = 1.25, 1.5, and 2
elseif Re_g >= 1E4
c_0 = [.245, .203, .162  ];   %b = 1.25, 1.5, and 2
c_1 = [.339E4, .248E4, .181E4 ]; %b = 1.25, 1.5, and 2
c_2 = [-.984E7, -.758E7, .792E8]; %b = 1.25, 1.5, and 2
c_3 = [.132E11, .104E11, -.165E13]; %b = 1.25, 1.5, and 2
c_4 = [-.599E13, -.482E13, .872E16];     %b = 1.25, 1.5, and 2
end 
% Assigning correct values to each tube spacing
% if (SL_sample == 1.25)
%     i = 1;
%     else if (SL_sample == 1.5)
%         i = 2;
%         else if (SL_sample == 2)
%             i = 3;
%             end
%         end
% end 

%%%%%%%%% Correction Factors %%%%%%%%%% 

%k_1 is the influence of pitch ratio
        %Will add later if data can be found

%k_2 =  Influence of Temperature on Fluid properties. Neglected 

% k_3 is Entry length effects 

%Entry loss coefficients (Re < 1E2 but < 1E4)
el_1 = [1.4, 1.3, 1.2, 1.1, 1, 1, 1];
%Entry loss coefficients (Re > 1E4 but < 1E6) 
el_2 = [1.1, 1.05, 1, 1, 1, 1, 1];
%Entry loss coefficients (Re > 1E6)
el_3 = [.25, .45, .6, .65, .7, .75, .8];

k_3 = 1; %Setting for tubes > 7
if (N_L < 7 && N_L > 0 && Re_g < 1E4 && Re_g > 1E2)
        k_3 = el_1(N_L);
elseif (N_L < 7 && N_L > 0 && Re_g < 1E6 && Re_g >= 1E4)
      k_3 = el_2(N_L);
elseif (N_L < 7 && N_L > 0 && Re_g > 1E6)
      k_3 = el_3(N_L);
end

%This section calculates Euler number at sample ratios and then
%uses linear interpolation for the actual ratio. 
SL_sample=[1.25,1.5,2];
Eu_p_sample=zeros(size(SL_sample));
Eu_sample=zeros(size(SL_sample));

for i=1:size(SL_sample,2)
%Power series for Euler number per row. From same website as above.
Eu_p_sample(i) = (c_0(i)./Re_g.^0)+(c_1(i)./Re_g.^1)+(c_2(i)./Re_g.^2)+(c_3(i)./Re_g.^3)+(c_4(i)./Re_g.^4);

%%%Corrected Euler Number    
Eu_sample(i) = Eu_p_sample(i).*k_3; %.*k_1;

end

Eu=interp1(SL_sample,Eu_sample,SL);

%Using the relation Eu = dP/((1/2)*rho*v^2)
dP = Eu.*((rho_g.*u_max_app.^2)./2); %Pressure drop per row 
dP_total = dP.*N_L/100000; %pressure drop across 10 rows expressed in bar
end
