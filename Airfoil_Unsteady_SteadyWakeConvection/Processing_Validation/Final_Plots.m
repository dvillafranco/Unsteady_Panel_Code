% This script will plot the necessary figures for paper

% Load necessary data
%load Run_Data/Unsteady_NACA0001_3000resolution.mat
%load Run_Data/Unsteady_NACA2412_3000steps.mat
%load Run_Data/NACA2412_FINAL_DATA.mat
%load  Run_Data/Unsteady_NACA0012.mat
load  Run_Data/NACA0012_good_data_final.mat

% Path of Vortex and Wake 
figure(1)
hp = plot(xnew,ynew,'k');
hold on
plot(x_imp,y_imp,'k',x_fixed,y_fixed,'k','LineWidth',1);
axis equal
title('Vortex and Wake Paths')


ts = 18;
fn = 'helvetica';
lw = 2;

set(gca, 'fontsize', ts)
set(gca, 'linewidth', lw)
set(hp, 'linewidth', lw)

set(gca, 'fontname', fn)
set(get(gca, 'XLabel'), 'FontSize', ts, 'FontName', fn)
set(get(gca, 'YLabel'), 'FontSize', ts, 'FontName', fn)
set(get(gca, 'ZLabel'), 'FontSize', ts, 'FontName', fn)
set(get(gca, 'title'), 'FontSize', ts, 'FontName', fn)
set(legend, 'Fontsize', ts, 'FontName', fn, 'linewidth', lw)

set(gca, 'gridlinestyle', '-')


figure(2)

% Establish Parameters 
gamma = 0.01; % Strength of Vortex
c = 1; %chord length
U = 1; % Freestream
rho = 1; % Density of air
h = 0.01; % Closest approach distance of vortex

F = ((real(fer2).^2+imag(fer2).^2).^0.5).*exp(-k*h)*pi*2;

hp = semilogy(pi*f,F);
xlim([0 25])
%ylim([10^-5 1]) 
xlabel('Normalized Frequency')
ylabel('Normalized Lift')

ts = 18;
fn = 'helvetica';
lw = 2;

set(gca, 'fontsize', ts)
set(gca, 'linewidth', lw)
set(hp, 'linewidth', lw)

set(gca, 'fontname', fn)
set(get(gca, 'XLabel'), 'FontSize', ts, 'FontName', fn)
set(get(gca, 'YLabel'), 'FontSize', ts, 'FontName', fn)
set(get(gca, 'ZLabel'), 'FontSize', ts, 'FontName', fn)
set(get(gca, 'title'), 'FontSize', ts, 'FontName', fn)
set(legend, 'Fontsize', ts, 'FontName', fn, 'linewidth', lw)

set(gca, 'gridlinestyle', '-')

figure(3)
hp = plot(time,lift_coef);
title('Lift')
xlim([0 12])
xlabel('Time (s)')
ylabel('Normalized Lift')

ts = 18;
fn = 'helvetica';
lw = 2;

set(gca, 'fontsize', ts)
set(gca, 'linewidth', lw)
set(hp, 'linewidth', lw)

set(gca, 'fontname', fn)
set(get(gca, 'XLabel'), 'FontSize', ts, 'FontName', fn)
set(get(gca, 'YLabel'), 'FontSize', ts, 'FontName', fn)
set(get(gca, 'ZLabel'), 'FontSize', ts, 'FontName', fn)
set(get(gca, 'title'), 'FontSize', ts, 'FontName', fn)
set(legend, 'Fontsize', ts, 'FontName', fn, 'linewidth', lw)

set(gca, 'gridlinestyle', '-')


 
