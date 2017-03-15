%%%%%%%%%%%%%%%%%%%%%%%%%% Panel_BEM_sears.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
% Author: Dorien O. Villafranco
% Unsteady Fluid Mechanics and Acoustics Lab
% Department of Mechanical Engineering
% College of Engineering
% Boston University
%
%
% This script will load the data from the unsteady panel code cascade case,
% the unsteady panel code single airfoil case and BEM single airfoil
% code. It will develop the Sears functions along with the analytical
% function for lift due to a passing vortex. The unsteady panel code data
% and the BEM data will be processed. This processing includes removing any
% irregularities in the data which are attributed to the trailing edge
% discretization of the airfoil. The FFT of both data sets will then be
% taken to obtain spectral results for the data sets. From this, plots of
% magnitude of lift vs. reduced frequency can be generated. The plots will
% show both the time signal from the BEM and panel code compared against
% each other along with the lift magnitude vs. reduced frequency for the
% panel code signal, BEM signal and the Sears lift signal. 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Clear and close all figure, set necessary parameters
clear
close all
rho = 1;

%% Load the Panel Code Data
%load L_unsteady_9x_heavymesh.mat
load cascade_h10.mat
c = 1;
%% Load the BEM Code data
D1 = BEM2d_readbin('2NACA0001_long');

%% Develop Sears Function & Lift
% Develop the Sears function 
k_1 = 0.01:0.01:100;
H0 = besselh(0,2,k_1);
H1 = besselh(1,2,k_1);
S = 2./(pi*k_1.*(H0-(1i*(H1))));
h = 0.02;
% Develop Lift scaled by \Gamma
L_2 = rho*(c/2).*(-1i*S).*exp(-k_1*h);
L_mag = abs(L_2);



%% Process the Panel Code Data (Interpolation & FFT)
% Calculate Phi Dot
h = 10;
phi_dot_casc;
time = 0.01:0.01:40;
% Interpolate the Signal 
%Lift_Panel_Interp;
Lift_Panel_Interp_4000steps;

% Conduct FFT and store panel FFT results
time_panel = 0.01:0.01:40;
%time_panel = 0.005:0.005:20;
foft = L_unsteady_fix;
ffer3;
close figure 1
fer_Panel = fer2;
f_Panel = f;

%% Process the BEM Code Data (Interpolation & FFT)
before_dip = 875;   
after_dip = 994;

izeropad =1.0;

% Set up the points for interpolation algorithm
L_front = D1.L(1:before_dip);
L_back = D1.L(after_dip:end);
time_good = [D1.t(1:before_dip) ; D1.t(after_dip:end)];
lift_good = [L_front, L_back];
time_query = D1.t(before_dip+1:after_dip-1);

% Interpolate points
L_int = interp1(time_good, lift_good, time_query,'pchip');

D1.L = [L_front, L_int',L_back];

% Extrapolate at the beginning
time_good = D1.t(6:end);
lift_good = D1.L(6:end);
time_query = D1.t(1:5);

L_ext = interp1(time_good,lift_good,time_query,'spline','extrap');
D1.L = [L_ext', lift_good];

% Conduct FFT and store BEM FFT results
clear time
time = D1.t;
foft = D1.L;
ffer3;
close figure 1
fer_BEM = fer2;
f_BEM = f;


%% Plot the Data 

% Plot the lift curve
figure(10)
hp1 = plot(time_panel,L_unsteady_fix);
hold on;
hp2 = plot(time,D1.L);
set(hp1,'LineWidth',2);
set(hp2,'LineWidth',2);
xlabel('Time (s)');
ylabel('Lift Coefficient');
grid on
set(gca,'FontSize',16);
leg1 = legend('Panel Code Signal','BEM Code Signal');
set(leg1,'Location','SouthEast');

% Plot the exact output of FFT
figure(20)
hp1 = semilogy(f_Panel,abs(fer_Panel));
hold on;
hp2 = semilogy(f_BEM,abs(fer_BEM));
leg2 = legend('Panel FFT Result','BEM FFT Result');
set(gca,'FontSize',16)
grid on;
xlabel('Frequency - f')
set(hp1,'LineWidth',1.5)
set(hp2,'LineWidth',1.5)


% Plot the scaled lift values vs. reduced frequency
figure(30)
hp1 = semilogy(f_Panel*pi,abs(fer_Panel/vortstrength)./exp(-f_Panel*h*pi/1)*2*pi,'o');
hold on;
hp2 = semilogy(f_BEM*pi,abs(fer_BEM/vortstrength)./exp(-f_BEM*h*pi)*pi,'o');
hp3 = semilogy(k_1,L_mag*2);
set(hp1,'LineWidth',1.5);
set(hp2,'LineWidth',1.5);
set(hp3,'LineWidth',1.5);
set(gca,'FontSize',15);
grid on;
xlabel('Reduced Frequency - k');
ylabel('Lift Coefficient');
xlim([ 0 20]);
leg3 = legend('Panel','BEM','Sears');


%% Ventres Comparison

% Set the Mach number
mach = 0.01;
h = 0.02;
% Import relevant data
filename = 'ventres_ksweep_h10.txt';
A = importdata(filename);

for i = 1:length(A)
    k_ven(i) = A(i,1);
    cl_ven(i) = A(i,2);
   
end
cl_ven = cl_ven.*exp(-k_ven*h);
figure(202)
hp = semilogy(k_ven,cl_ven,'--');
hold on;
hp1 = semilogy(f_Panel*pi,abs(fer_Panel/vortstrength*2*pi))%.*exp(-f_Panel*h*pi/1)*2*pi);