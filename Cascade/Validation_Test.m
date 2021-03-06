%% Validation Method 
% as Presented in Unsteady Blade Response:The BVI Model vs. 
% the Gust Model (Grace, 2000)

% Analytical Lift Response of a Flat Plate to a passing vortex

% Establish normalized frequency k = (wc/2U) vector
k = linspace(0,1000,84);


% Establish Parameters 
gamma = 0.01; % Strength of Vortex
c = 1; %chord length
U = 1; % Freestream
rho = 1; % Density of air
h = 0.01; % Closest approach distance of vortex

F = ((real(fer2).^2+imag(fer2).^2).^0.5).*exp(-k*0.02)*2*pi;

H0 = besselh(0,2,k);
H1 = besselh(1,2,k);
S = 2./(pi*k.*(H0-(1i*(H1))));
L = rho*gamma*(c/2).*exp((-1*k*h)).*(-1i*S)*(2*pi);
L_mag = (real(L).^2 + imag(L).^2).^0.5;

% figure(100)
% plot(k/pi,L_mag,(f),(real(fer2).^2+imag(fer2).^2).^0.5,'x','LineWidth',1)
% ylabel('Normalized Lift')
% xlabel('Normalized Frequency')
% title('Lift Response of Flat Plate to Passing Vortex')
% legend('Analytical','BVI Data')
% xlim([0 15])
% %set(gca,'fontsize',14)
% grid on

figure(23)
semilogy(k/pi,L_mag,f,F,'o')
xlim([0 300])
ylabel('Normalized Lift')
xlabel('Normalized Frequency')
title('Lift Response of Cascade of NACA 0012 to Passing Vortex')
legend('Analytical','BVI Data')
grid on



% 
% %% Validation Method 
% % as Presented in Unsteady Blade Response:The BVI Model vs. 
% % the Gust Model (Grace, 2000)
% 
% % Analytical Lift Response of a Flat Plate to a passing vortex
% 
% % Establish normalized frequency k = (wc/2U) vector;
% k = linspace(0,500,85);
% 
% 
% % Establish Parameters 
% gamma = 0.01; % Strength of Vortex
% c = 1; %chord length
% % sos = 340;
% % mach = 0.2;
% % U = mach*sos;
% U = 1; %Freestream
% rho = 1.225; % Density of air
% h = 0.08; % Closest approach distance of vortex
% 
% F = ((real(fer2/gamma).^2+imag(fer2/gamma).^2).^0.5).*exp(-k*h)/2/pi;%*(2*pi)^4;
% 
% H0 = besselh(0,2,k);
% H1 = besselh(1,2,k);
% S = 2./(pi*k.*(H0-(1i*(H1))));
% %S = ((2*i)./(pi*k)).*(exp(-i*k)./(H1+i*H0));
% L = rho*gamma*(c/2).*exp((-1*k*h)).*(-1i*S)*2*pi;%.*exp(k*h);
% L_mag = (real(L).^2 + imag(L).^2).^0.5;
% 
% % figure(100)
% % plot(k/pi,L_mag,(f),(real(fer2).^2+imag(fer2).^2).^0.5,'x','LineWidth',1)
% % ylabel('Normalized Lift')
% % xlabel('Normalized Frequency')
% % title('Lift Response of Flat Plate to Passing Vortex')
% % legend('Analytical','BVI Data')
% % xlim([0 15])
% % %set(gca,'fontsize',14)
% % grid on
% 
% figure(23)
% semilogy(k/pi,L_mag,f,F)
% xlim([0 100])
% ylabel('Normalized Lift')
% xlabel('Normalized Frequency')
% title('Lift Response of Flat Plate to Passing Vortex')
% legend('Analytical','BVI Data')
% grid on
% set(gca,'fontsize',14)
% 
% 
% 
% %%%%%%% Validation for Unsteady Thickness Case - Incompressible
% t_max = 0.06;
% S = (1./(1+(2*pi*k).^1.3)).^(1/1.3);
% x_c = [0.01,0.1];
% tau_c = zeros(1,2);
% for i= 1:2;
% tau_c(i) = t_max/0.2*(0.2969*(x_c(i)^0.5)-0.1260*(x_c(i)^2)+0.2843*(x_c(i)^3)-...
%     0.1015*(x_c(i)^4));
% end
% 
% D_sum = sum(tau_c);
% 
% beta_atten = (2*pi*log(10)/10)*(1.123*D_sum+5.317*D_sum^2);
% 
% H_resp = (abs(S)).^2.*exp(-1*beta_atten*k/pi)./(pi*rho*c*1)^2*(1+0.8*t_max)^2;
% H_mag = (real(H_resp).^2+imag(H_resp).^2).^0.5;
% 
% figure(33)
% semilogy(k/pi,H_mag,f,F,'*')
% xlim([0 100])
% grid on
% ylabel('Normalized Lift')
% xlabel('Normalized Frequency')
% set(gca,'FontSize',14)
% 
% 
% 
