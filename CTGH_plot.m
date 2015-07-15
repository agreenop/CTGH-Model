%This function will map the temperatures, pressures, and heat transfer
%outputs for the liquid coolant and gas medium for the CTGH.
function CTGH_plot(T_l,T_g,Q,P_l,P_g)
%Plot liquid coolant temperatures
rho=10:size(T_l,1)+9;
theta = linspace(90,-270,size(T_l,2))*pi/180;
[th, r] = meshgrid(theta, rho);
figure(1)
surf(r.*cos(th),r.*sin(th),T_l,'linestyle','none'); 
view(2);
axis equal tight;
colorbar;
%Plot gas temperatures
rho_g=10:size(T_g,1)+9;
theta_g = linspace(90,-270,size(T_g,2))*pi/180;
[th_g, r_g] = meshgrid(theta_g, rho_g);
figure(2)
surf(r_g.*cos(th_g),r_g.*sin(th_g),T_g,'linestyle', 'none'); 
view(2);
axis equal tight;
colorbar;
%Plot heat transfer
figure(3)
surf(r.*cos(th),r.*sin(th),Q,'linestyle', 'none'); 
view(2);
axis equal tight;
colorbar;
%Plot coolant pressure drop
figure(4)
surf(r.*cos(th),r.*sin(th),P_l,'linestyle', 'none'); 
view(2);
axis equal tight;
colorbar;
%Plot gas pressure drop
figure(5)
surf(r_g.*cos(th_g),r_g.*sin(th_g),P_g,'linestyle', 'none'); 
view(2);
axis equal tight;
colorbar;