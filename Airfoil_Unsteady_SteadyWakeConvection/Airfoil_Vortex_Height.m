

clear 
close all

time = 0.01:0.01:40;

% Load Airfoil vortex height h = 0.2
load Panel_Unsteady_hvor02.mat
phi_dot_calc;
L_unsteady_hvor02 = L_unsteady;
h = 0.2;
%Lift_Panel_Interp_4000steps;
%L_unsteady02 = L_unsteady_fix;
% foft = L_unsteady02;
foft = L_unsteady_hvor02;
ffer3;
close figure 1
fer_Panel = fer2;
f_Panel = f;
figure(11111)
hp2 = semilogy(f_Panel*pi,abs(fer_Panel/vortstrength));%./exp(-f_Panel*h*pi/1)*2*pi,'o');
hold on;


% Load Airfoil vortex height h = 0.3
clear L_unsteady foft
load Panel_Unsteady_hvor03.mat
phi_dot_calc;
L_unsteady_hvor03 = L_unsteady;
h = 0.3;
%Lift_Panel_Interp_4000steps;
%L_unsteady02 = L_unsteady_fix;
%foft = L_unsteady02;
foft = L_unsteady_hvor03;
ffer3;
close figure 1
fer_Panel = fer2;
f_Panel = f;
figure(11111)
hp3 = semilogy(f_Panel*pi,abs(fer_Panel/vortstrength));%./exp(-f_Panel*h*pi/1)*2*pi,'o');
hold on;


% Load airfoil vortex height h = 0.4


% Load airfoil vortex height h = 0.5
clear L_unsteady;
load Panel_Unsteady_hvor05.mat
phi_dot_calc;
L_unsteady_hvor05 = L_unsteady;
h = 0.5;
%Lift_Panel_Interp_4000steps;
%L_unsteady02 = L_unsteady_fix;
%foft = L_unsteady02;
foft = L_unsteady_hvor05;
ffer3;
close figure 1
fer_Panel = fer2;
f_Panel = f;
figure(11111)
hp5 = semilogy(f_Panel*pi,abs(fer_Panel/vortstrength));%./exp(-f_Panel*h*pi/1)*2*pi,'o');
hold on;


% % Load airfoil vortex height h = 0.6
% clear L_unsteady;
% load Panel_Unsteady_hvor06.mat
% phi_dot_calc;
% L_unsteady_hvor06 = L_unsteady;


% Load airfoil vortex height h = 0.7
clear L_unsteady;
load Panel_Unsteady_hvor07.mat
phi_dot_calc;
L_unsteady_hvor07 = L_unsteady;
h = 0.7;
%Lift_Panel_Interp_4000steps;
%L_unsteady02 = L_unsteady_fix;
%foft = L_unsteady02;
foft = L_unsteady_hvor07;
ffer3;
close figure 1
fer_Panel = fer2;
f_Panel = f;
figure(11111)
hp7 = semilogy(f_Panel*pi,abs(fer_Panel/vortstrength));%./exp(-f_Panel*h*pi/1)*2*pi,'o');
hold on;


% Load airfoil vortex height h = 0.8
clear L_unsteady;
load Panel_Unsteady_hvor08.mat;
phi_dot_calc;
L_unsteady_hvor08 = L_unsteady;
h = 0.8;
foft = L_unsteady_hvor08;
ffer3;
close figure 1
fer_Panel = fer2;
f_Panel = f;
figure(11111)
hp8 = semilogy(f_Panel*pi,abs(fer_Panel/vortstrength));%./exp(-f_Panel*h*pi/1)*2*pi,'o');
hold on;

% % Load airfoil vortex height h = 0.9
% clear L_unsteady
% load Panel_Unsteady_hvor09.mat;
% phi_dot_calc;
% L_unsteady_hvor09 = L_unsteady;

% Load airfoil vortex height h = 0.1
clear L_unsteady
load Panel_Unsteady_hvor1.mat
phi_dot_calc;
L_unsteady_hvor1 = L_unsteady;
h = 0.9;
% Lift_Panel_Interp_4000steps;
% L_unsteady02 = L_unsteady_fix;
foft = L_unsteady_hvor1;
ffer3;
close figure 1
fer_Panel = fer2;
f_Panel = f;
figure(11111)
hp10 = semilogy(f_Panel*pi,abs(fer_Panel/vortstrength));%./exp(-f_Panel*h*pi/1)*2*pi,'o');
xlim([0 20]);
grid on;
set(gca,'FontSize',15);
legend('h = 0.2','h = 0.3','h = 0.5','h = 0.7','h = 0.8', 'h = 1.0');
set(hp2,'LineWidth',2);
set(hp3,'LineWidth',2);
%set(hp4,'LineWidth',2);
set(hp5,'LineWidth',2);
% set(hp6,'LineWidth',2);
set(hp7,'LineWidth',2);
set(hp8,'LineWidth',2);
%set(hp9,'LineWidth',2);
set(hp10,'LineWidth',2);

% Plot the various lift signals. 
figure(1)
hp2 = plot(time,L_unsteady_hvor02);
hold on;
hp3 = plot(time,L_unsteady_hvor03);
%hp4 = plot(time,L_unsteady_hvor04);
hp5 = plot(time,L_unsteady_hvor05);
%hp6 = plot(time,L_unsteady_hvor06);
hp7 = plot(time,L_unsteady_hvor07);
hp8 = plot(time,L_unsteady_hvor08);
% hp9 = plot(time,L_unsteady_hvor09);
hp10 = plot(time,L_unsteady_hvor1);
set(gca,'FontSize',15);
set(hp2,'LineWidth',2)
set(hp3,'LineWidth',2)
%set(hp4,'LineWidth',2)
set(hp5,'LineWidth',2)
%set(hp6,'LineWidth',2)
set(hp7,'LineWidth',2)
set(hp8,'LineWidth',2)
% set(hp9,'LineWidth',2)
set(hp10,'LineWidth',2)
xlabel('Time (s)')
ylabel('Lift Coefficient')
legend('h = 0.2','h = 0.3','h = 0.5','h = 0.7','h = 0.8', 'h = 1.0');
xlim([0.5 11]);