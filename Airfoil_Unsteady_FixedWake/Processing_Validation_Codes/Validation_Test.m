%% Validation Method 
% as Presented in Unsteady Blade Response:The BVI Model vs. 
% the Gust Model (Grace, 2000)

% Analytical Lift Response of a Flat Plate to a passing vortex

% Establish normalized frequency k = (wc/2U) vector
k = linspace(0,50,600);


% Establish Parameters 
gamma = 0.02; % Strength of Vortex
c = 1; %chord length
U = 1; % Freestream
rho = 1; % Density of air
h = 0.1; % Closest approach distance of vortex

F = ((real(fer2).^2+imag(fer2).^2).^0.5).*exp(k*0.1);

H0 = besselh(0,2,k);
H1 = besselh(1,2,k);
S = 2./(pi*k.*(H0-(1i*(H1))));
L = rho*gamma*(c/2).*exp((-1*k*h)).*(-1i*S)*(2*pi);
L_mag = (real(L).^2 + imag(L).^2).^0.5;

figure(100)
plot(k/pi,L_mag,(f),(real(fer2).^2+imag(fer2).^2).^0.5,'x','LineWidth',1)
ylabel('Normalized Lift')
xlabel('Normalized Frequency')
title('Lift Response of Flat Plate to Passing Vortex')
legend('Analytical','BVI Data')
xlim([0 15])
%set(gca,'fontsize',14)
grid on

figure(23)
semilogy(k/pi,L_mag,f,F)
xlim([0 15])
ylabel('Normalized Lift')
xlabel('Normalized Frequency')
title('Lift Response of Flat Plate to Passing Vortex')
legend('Analytical','BVI Data')
grid on