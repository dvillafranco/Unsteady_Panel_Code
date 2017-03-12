% This script will show the difference between the two results for the
% Panel code where there is convection in one due to to the local
% velocities of the vorticies in the wake, and the imposed vortex. For the
% other data set, all the velocites were calculated in the steady case. In
% this case, on the freestream velocity anc velocity contribution due to
% the linearly varying panels were considered. 


% Load the Fast data (steady velocity convection)
Panel_BEM_Sears;



% Load the Slow data (unsteady velocity convection)
clear; 
close figure 20

load('L_unsteady_9x_heavymesh.mat')
phi_dot_calc;
time = 0.01:0.01:20;
Lift_Panel_Interp;
foft = L_unsteady_fix;
ffer3;
fer_unsteady_conv = fer2;

% Plot all the lift signals 
figure(10)
plot(time,L_unsteady_fix,'LineWidth',2);
legend('Panel - steady convection','BEM','Panel - unsteady convection');

figure(30);
hp11 = semilogy(f*pi,abs(fer_unsteady_conv/vortstrength)./exp(-f*h*pi/1)*pi);