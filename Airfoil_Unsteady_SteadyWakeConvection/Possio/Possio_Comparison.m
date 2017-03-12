% This code will create the figure which compares the BVI calculation for a
% thin airfoil ~ NACA 0002 to the Possio calculation for a Mach number of
% 0.5. 
k = linspace(0,60,600);
h = 0.18;

%load Unsteady_Possio_Mach03.mat

load Unsteady_Comp_Possio_mach05.mat

figure(1)
hp = semilogy(k,lift.*exp(-k*h)/4/pi^2);
xlim([0 20])

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



load Unsteady_Good_BVIMach03

h = 0.08;
% Establish normalized frequency k = (wc/2U) vector
k = linspace(0,60,600);


% Establish Parameters 
gamma = 0.02; % Strength of Vortex
c = 1; %chord length
sos = 340;
mach = 0.5;
U = mach*sos;
%U = 1; %Freestream
rho = 1.0; % Density of air
 % Closest approach distance of vortex

F = ((real(fer2/gamma).^2+imag(fer2/gamma).^2).^0.5).*exp(-k*0.01)*100/2/pi;
hold on
hp = semilogy(k,F,'--','LineWidth',2);

legend('Possio','BVI')
xlabel('Normalized Frequency')
ylabel('Normalized Lift')

