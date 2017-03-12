% This script will take the unsteady lift signal from the panel code and
% perform the interpolation necessary to remove the unphysical spike at the
% trailing edge
time = 0.005:0.005:20;
rho = 1;
c = 1;
h = 0.02;

% Find the points to the left of the spike which are to be kept, and the
% points to the right which are to be kept.
before_dip      = 1728;    %400; 500    906
after_dip       = 1940;   %530;  620   1061
 
% % Set up the points for interpolation algorithm
L_front = L_unsteady(1:before_dip);
L_back  = L_unsteady(after_dip:end);
time_good = [time(1:before_dip) , time(after_dip:end)];
lift_good = [L_front, L_back]; 
time_query = time(before_dip+1:after_dip-1);

% % Interpolate points
L_int = interp1(time_good, lift_good, time_query,'pchip');
L_unsteady_fix = [L_front, L_int,L_back];
% 
% before_dip = 1036
% after_dip = 1099
% 
% % Set up the points for interpolation algorithm
% L_front = L_unsteady_fix(1:before_dip);
% L_back  = L_unsteady_fix(after_dip:end);
% time_good = [time(1:before_dip) , time(after_dip:end)];
% lift_good = [L_front, L_back]; 
% time_query = time(before_dip+1:after_dip-1);
% 
% % Interpolate points
% L_int = interp1(time_good, lift_good, time_query,'pchip');
% L_unsteady_fix = [L_front, L_int,L_back];
% 
% % Interpolate 3
% before_dip = 1339
% after_dip = 1365
% 
% % Set up the points for interpolation algorithm
% L_front = L_unsteady_fix(1:before_dip);
% L_back  = L_unsteady_fix(after_dip:end);
% time_good = [time(1:before_dip) , time(after_dip:end)];
% lift_good = [L_front, L_back]; 
% time_query = time(before_dip+1:after_dip-1);
% 
% % Interpolate points
% L_int = interp1(time_good, lift_good, time_query,'pchip');
% L_unsteady_fix = [L_front, L_int,L_back];

% Extrapolate at the beginning
time_good = time(31:end);
lift_good = L_unsteady_fix(31:end);
time_query = time(1:30);

L_ext = interp1(time_good,lift_good,time_query,'pchip','extrap');
L_unsteady_fix = [L_ext, lift_good];

% %% Plotting
% 
% figure(2)
% hp3 = plot(time,L_unsteady);
% hold on;
% hp4 = plot(time,L_unsteady_fix);
% set(gca,'FontSize',15);
% grid on;
% 
% hp5 = plot(D1.t,D1.L);
% set(hp3,'LineWidth',2);
% set(hp4,'LineWidth',2);
% set(hp5,'LineWidth',2);
% xlabel('Time (s)');
% ylabel('Lift Coefficient');
% leg2 = legend('Original Signal','Interpolated Signal','BEM Signal');
% set(leg2,'FontSize',13);
% xlim([0 10])

% 
% 
% figure(1000)
% hp1 = semilogy(f*pi,abs(fer2/vortstrength)./exp(-f*h*pi/1)*pi);
% %hp1 = semilogy(f*pi,abs(fer2/vortstrength)*pi);
% hold on
% hps = semilogy(k_1,abs(L_mag)*2);
% 
% foft = D1.L;
% time = 0.0125:0.0125:10;
% ffer3;
% close figure 1
% figure(1000)
% hp2 = semilogy(f*pi,abs(fer2/vortstrength)./exp(-f*h*pi/1)*pi/2);
% 
% 
