

clear 
close all

% Load Vortext height h = 0.02
load Cascade_h1_hvor02.mat

L_unsteady02_orig = L_unsteady;
h = 0.02;
Lift_Panel_Interp_4000steps;
L_unsteady02 = L_unsteady_fix;
foft = L_unsteady02;
ffer3;
close figure 1
fer_Panel = fer2;
f_Panel = f;
figure(11111)
hp1 = semilogy(f_Panel*pi,abs(fer_Panel/vortstrength),'>');%./exp(-f_Panel*h*pi/1)*2*pi,'o');
hold on;

% % Load Vortex height h = 0.04;
clearvars -except L_unsteady02_orig L_unsteady02
load cascade_fp_h1_h04.mat
L_unsteady04_orig = L_unsteady;
h = 0.04;
Lift_Panel_Interp_4000steps;
L_unsteady04 = L_unsteady_fix;
foft = L_unsteady_fix;
ffer3;
close figure 1
fer_Panel = fer2;
f_Panel = f;
figure(11111)
hp2 = semilogy(f_Panel*pi,abs(fer_Panel/vortstrength))%./exp(-f_Panel*h*pi/1)*2*pi);
% 
% % Load Vortex height h = 0.06
clearvars -except L_unsteady02_orig L_unsteady02 L_unsteady04_orig L_unsteady04
load cascade_fp_h1_h06.mat
phi_dot_casc;
L_unsteady06_orig = L_unsteady;

h = 0.06;
Lift_Panel_Interp_4000steps;
L_unsteady06 = L_unsteady_fix;
foft = L_unsteady_fix;
ffer3;
close figure 1
fer_Panel = fer2;
f_Panel = f;
figure(11111)
hp3 = semilogy(f_Panel*pi,abs(fer_Panel/vortstrength))%./exp(-f_Panel*h*pi/1)*2*pi);
% 
% Load Vortex Height h = 0.08
clearvars -except L_unsteady02_orig L_unsteady02 L_unsteady04_orig L_unsteady04 L_unsteady06_orig L_unsteady06
load cascade_fp_h1_h08.mat
L_unsteady08_orig = L_unsteady;
h = 0.08;
Lift_Panel_Interp_4000steps;
L_unsteady08 = L_unsteady_fix;
foft = L_unsteady_fix;
ffer3;
close figure 1
fer_Panel = fer2;
f_Panel = f;
figure(11111)
hp4 = semilogy(f_Panel*pi,abs(fer_Panel/vortstrength))%./exp(-f_Panel*h*pi/1)*2*pi,'k');

% Load vortex height h = 0.1
clearvars -except L_unsteady02_orig L_unsteady02 L_unsteady04_orig L_unsteady04 L_unsteady06_orig L_unsteady06 L_unsteady08_orig L_unsteady08
load cascade_fp_h1_h1.mat
L_unsteady1_orig = L_unsteady;
h = 0.1;
Lift_Panel_Interp_4000steps;
L_unsteady1 = L_unsteady_fix;
foft = L_unsteady_fix;
ffer3;
close figure 1
fer_Panel = fer2;
f_Panel = f;
figure(11111)
hp5 = semilogy(f_Panel*pi,abs(fer_Panel/vortstrength))%./exp(-f_Panel*h*pi/1)*2*pi);

%Load vortex height h = 0.2
clear L_unsteady_fix;
load cascade_fp_h1_h2.mat
L_unsteady2_orig = L_unsteady;
h = 0.2;
Lift_Panel_Interp_4000steps;
L_unsteady1 = L_unsteady_fix;
foft = L_unsteady_fix;
ffer3;
close figure 1
fer_Panel = fer2;
f_Panel = f;
figure(11111)
hp5 = semilogy(f_Panel*pi,abs(fer_Panel/vortstrength))%./exp(-f_Panel*h*pi/1)*2*pi,'o');

%Load vortex height h = 0.3
clear L_unsteady_fix;
load cascade_fp_h1_h3.mat
L_unsteady3_orig = L_unsteady;
h = 0.3;
Lift_Panel_Interp_4000steps;
L_unsteady1 = L_unsteady_fix;
foft = L_unsteady_fix;
ffer3;
close figure 1
fer_Panel = fer2;
f_Panel = f;
figure(11111)
hp6 = semilogy(f_Panel*pi,abs(fer_Panel/vortstrength))%./exp(-f_Panel*h*pi/1)*2*pi,'o');

%Load vortex height h = 0.4
clear L_unsteady_fix;
load cascade_fp_h1_h4.mat
L_unsteady4_orig = L_unsteady;
h = 0.4;
Lift_Panel_Interp_4000steps;
L_unsteady1 = L_unsteady_fix;
foft = L_unsteady_fix;
ffer3;
close figure 1
fer_Panel = fer2;
f_Panel = f;
figure(11111)
hp7 = semilogy(f_Panel*pi,abs(fer_Panel/vortstrength))%./exp(-f_Panel*h*pi/1)*2*pi,'o');

%Load vortex height h = 0.5
clear L_unsteady_fix;
load cascade_fp_h1_h5.mat
L_unsteady5_orig = L_unsteady;
h = 0.5;
Lift_Panel_Interp_4000steps;
L_unsteady1 = L_unsteady_fix;
foft = L_unsteady_fix;
ffer3;
close figure 1
fer_Panel = fer2;
f_Panel = f;
figure(11111)
hp8 = semilogy(f_Panel*pi,abs(fer_Panel/vortstrength));%./exp(-f_Panel*h*pi/1)*2*pi,'o');

%Load vortex height h = 0.6
clear L_unsteady_fix;
load cascade_fp_h1_h6.mat
L_unsteady6_orig = L_unsteady;
h = 0.6;
Lift_Panel_Interp_4000steps;
L_unsteady1 = L_unsteady_fix;
foft = L_unsteady_fix;
ffer3;
close figure 1
fer_Panel = fer2;
f_Panel = f;
figure(11111)
hp9 = semilogy(f_Panel*pi,abs(fer_Panel/vortstrength));%./exp(-f_Panel*h*pi/1)*2*pi,'o');

%Load vortex height h = 0.7
clear L_unsteady_fix;
load cascade_fp_h1_h4.mat
L_unsteady7_orig = L_unsteady;
h = 0.7;
Lift_Panel_Interp_4000steps;
L_unsteady1 = L_unsteady_fix;
foft = L_unsteady_fix;
ffer3;
close figure 1
fer_Panel = fer2;
f_Panel = f;
figure(11111)
hp10 = semilogy(f_Panel*pi,abs(fer_Panel/vortstrength));%./exp(-f_Panel*h*pi/1)*2*pi,'o');

%Load vortex height h = 0.8
clear L_unsteady_fix;
load cascade_fp_h1_h8.mat
L_unsteady8_orig = L_unsteady;
h = 0.8;
Lift_Panel_Interp_4000steps;
L_unsteady1 = L_unsteady_fix;
foft = L_unsteady_fix;
ffer3;
close figure 1
fer_Panel = fer2;
f_Panel = f;
figure(11111)
hp11 = semilogy(f_Panel*pi,abs(fer_Panel/vortstrength));%./exp(-f_Panel*h*pi/1)*2*pi,'o');

%Load vortex height h = 0.9
clear L_unsteady_fix;
load cascade_fp_h1_h9.mat
L_unsteady9_orig = L_unsteady;
h = 0.9;
Lift_Panel_Interp_4000steps;
L_unsteady1 = L_unsteady_fix;
foft = L_unsteady_fix;
ffer3;
close figure 1
fer_Panel = fer2;
f_Panel = f;
figure(11111)
hp12 = semilogy(f_Panel*pi,abs(fer_Panel/vortstrength),'d');%./exp(-f_Panel*h*pi/1)*2*pi,'o');


mach = 0.01;
h = 0.02;
% Import relevant data
filename = 'ventres_ksweep_h1.txt';
A = importdata(filename);

for i = 1:length(A)
    k_ven(i) = A(i,1);
    cl_ven(i) = A(i,2);
   
end
cl_ven = cl_ven%.*exp(-k_ven*h);
figure(11111)
hp = semilogy(k_ven,cl_ven);
grid on;
leg1 = legend('h = 0.02','h = 0.04','h = 0.06','h = 0.08','h = 0.1','h = 0.2','h=0.3','h=0.4','h=0.5','h=0.6','h=0.7','h=0.8','h=0.9','Ventres');
set(gca,'FontSize',15)
xlim([0 20])
xlabel('Reduced Frequency')
ylabel('FFT of lift coef')

figure(2)
plot(L_unsteady02);
hold on;
plot(L_unsteady04);
plot(L_unsteady06);
plot(L_unsteady08);
plot(L_unsteady1);

figure(3)
plot(L_unsteady02_orig);
hold on;
plot(L_unsteady04_orig);
plot(L_unsteady06_orig);
plot(L_unsteady08_orig);
plot(L_unsteady1_orig);


