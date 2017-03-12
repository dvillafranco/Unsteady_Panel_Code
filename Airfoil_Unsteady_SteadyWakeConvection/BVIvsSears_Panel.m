% This script will compare the BVI lift spectrum for a flat plate to the
% Sears function developed for a sinusoidal gust. 
clear all
close all
% Establish variables
c = 1;          % Chord length of airfoil
rho = 1;        % Density of flow
gamma = -0.1;    % Strength of imposed vortex
h = 0.1;
time = 0.0125:0.0125:10;
% load Unsteady_NACA0001_convect.mat
% c = 1;          % Chord length of airfoil

% Develop the Sears function 
k_1 = 0.01:0.01:100;
H0 = besselh(0,2,k_1);
H1 = besselh(1,2,k_1);
S = 2./(pi*k_1.*(H0-(1i*(H1))));
L_2 = rho*(c/2)*(-1i*S).*exp(-k_1*h);
L_mag = abs(L_2);


% Zero pad time signal to the left
% left_add = lift_coef(1);
% left_vec = left_add*ones(1,1000);
% left_vec = zeros(1,1000);
% lift_coef = [left_vec,lift_coef];
% %time = 0.0125:0.0125:10;
time = 0.0125:0.0125:10;
% time = linspace(0,8,length(lift_coef));
% 
% before_dip = 990;
% after_dip =  1300;
% 
% % Set up the points for interpolation algorithm
% L_front = lift_coef(1:before_dip);
% L_back = lift_coef(after_dip:end);
% time_good = [time(1:before_dip) , time(after_dip:end)];
% lift_good = [L_front, L_back]; 
% time_query = time(before_dip+1:after_dip-1);
% 
% % Interpolate points
% L_int = interp1(time_good, lift_good, time_query,'spline');
% lift_coef = [L_front, L_int,L_back];
% 
% 
% 
% % Zero pad time signal to the right
% right_add = lift_coef(end);
% right_vec = right_add*ones(1,1000);
% lift_coef = [lift_coef,right_vec];
% 
% % Interpolate the right side of the lift signal
% time = linspace(0,8,length(lift_coef));
% 
% before_dip = 1750;
% after_dip =  1850;
% 
% % Set up the points for interpolation algorithm
% L_front = lift_coef(1:before_dip);
% L_back = lift_coef(after_dip:end);
% time_good = [time(1:before_dip) , time(after_dip:end)];
% lift_good = [L_front, L_back]; 
% time_query = time(before_dip+1:after_dip-1);
% 
% % Interpolate points
% L_int = interp1(time_good, lift_good, time_query,'spline');
% lift_coef = [L_front, L_int,L_back];

% figure(2000)
% plot(time,lift_coef*2)
% hold on;
% D1 = BEM2d_readbin('2NACA0001_shift');
% plot(D1.t,D1.L,'--','LineWidth',2);
% legend('Panel Lift','BEM Lift')


% Develop the Sears function 
k_1 = 0.01:0.01:100;
H0 = besselh(0,2,k_1);
H1 = besselh(1,2,k_1);
S = 2./(pi*k_1.*(H0-(1i*(H1))));
L_2 = rho*(c/2)*(-1i*S).*exp(-k_1*h);
L_mag = abs(L_2);
foft = lift_coef;
ffer3;
close figure 1
figure(1000)
hp1 = semilogy(f*pi,abs(fer2/vortstrength)./exp(-f*h*pi/1)*pi,'x');
%hp1 = semilogy(f*pi,abs(fer2/vortstrength)*pi);
hold on
hps = semilogy(k_1,abs(L_mag)*2);
xlim([0 20])
foft = D1.L;
time = D1.t;
ffer3;
figure(1000)
hp24=semilogy(f*pi,abs(fer2/vortstrength)./exp(-f*h*pi/1)*pi/2,'--');


load Unsteady_NACA0001_convect.mat
time = 0.0125:0.0125:10;figure(2000)
plot(time,lift_coef*2,'x')

load Unsteady_NACA0001_convect.mat
figure(3000)
plot(xnew,ynew,'-');
hold on;
plot(D1.Xc(:,1),D1.Xc(:,2));
plot(x_imp,y_imp,'o','LineWidth',2);
Vort_2NACA0001
plot(Xvort(:,1),Xvort(:,2),'x');
axis equal
plot(x_fixed,y_fixed,'o','LineWidth',2);
plot(D1.Xw{1,800}(:,1),D1.Xw{1,800}(:,2),'x');
legend('Panel Body','BEM Body','Panel Vortex','BEM Vortex','Panel Wake','BEM Wake')
