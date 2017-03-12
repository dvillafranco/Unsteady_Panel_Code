% This script will show the plots for Sears, BVI, Thickness Correction and
% Thickness + Camber


% Load Necessary Data

load Run_Data/NACA0012_good_data_final.mat

%BVI Formulation 
% Establish Parameters 
gamma = 0.01; % Strength of Vortex
c = 1; %chord length
U = 1; % Freestream
rho = 1; % Density of air
h = 0.08; % Closest approach distance of vortex

F = ((real(fer2).^2+imag(fer2).^2).^0.5).*exp(-k*h);


% Sears Forumlation
H0 = besselh(0,2,k);
H1 = besselh(1,2,k);
S = 2./(pi*k.*(H0-(1i*(H1))));
%S = ((2*i)./(pi*k)).*(exp(-i*k)./(H1+i*H0));
L = rho*gamma*(c/2).*exp((-1*k*h)).*(-1i*S);%.*exp(-k*h);
L_mag = (real(L).^2 + imag(L).^2).^0.5;

% Thickness Formulation 
t_max = 0.12;
S = (1./(1+(2*pi*k).^1.3)).^(1/1.3);
x_c = [0.01,0.1];
tau_c = zeros(1,2);
for i= 1:2;
tau_c(i) = t_max/0.2*(0.2969*(x_c(i)^0.5)-0.1260*(x_c(i)^2)+0.2843*(x_c(i)^3)-...
    0.1015*(x_c(i)^4));
end

D_sum = sum(tau_c);

beta_atten = (2*pi*log(10)/10)*(1.123*D_sum+5.317*D_sum^2);

H_resp = (abs(S)).^2.*exp(-1*beta_atten*k/pi)./(pi*rho*c*1)^2*(1+0.8*t_max)^2;
H_mag = (real(H_resp).^2+imag(H_resp).^2).^0.5.*exp(-h*k);


figure(1)
hp2 =semilogy(k,L_mag);

hold on
hp = semilogy(pi*f,F,'--');
hp3 = semilogy(k,H_mag,'-.');
xlim([0 20])
%ylim([10^-5 1]) 
xlabel('Normalized Frequency')
ylabel('Normalized Lift')

%load Run_Data/Unsteady_NACA2412.mat
load Run_Data/NACA2412_FINAL_DATA.mat
%BVI Formulation 
% Establish Parameters 
gamma = 0.01; % Strength of Vortex
c = 1; %chord length
U = 1; % Freestream
rho = 1; % Density of air
h = 0.2; % Closest approach distance of vortex

F = ((real(fer2).^2+imag(fer2).^2).^0.5).*exp(-k*h);

hp4 = semilogy(k,F,':');

ts = 18;
fn = 'helvetica';
lw = 3;
lw2 = 2;
lw3 = 3;
lw4 = 3;

set(gca, 'fontsize', ts)
set(gca, 'linewidth', lw2)
set(hp, 'linewidth', lw)
set(hp2, 'linewidth',lw2)
set(hp3,'linewidth',lw3)
set(hp4,'linewidth',lw4)

set(gca, 'fontname', fn)
set(get(gca, 'XLabel'), 'FontSize', ts, 'FontName', fn)
set(get(gca, 'YLabel'), 'FontSize', ts, 'FontName', fn)
set(get(gca, 'ZLabel'), 'FontSize', ts, 'FontName', fn)
set(get(gca, 'title'), 'FontSize', ts, 'FontName', fn)
set(legend, 'Fontsize', ts, 'FontName', fn, 'linewidth', lw2)

set(gca, 'gridlinestyle', '-')
legend('Sears','BVI Thickness','Analytical Thickness','BVI Thickness & Camber')
