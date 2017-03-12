% Dorien Villafranco
% Unsteady Aerodynamics and Aeroacoustics Lab
% Department of Mechanical Engineering
% Boston University 
% Janurary 06, 2016

% This script will serve as validation for the 2D incompressible unsteady 
% flat plate vortex shedding problem 


% Define U Infinity
U = 1;

% Define chord length
c = 1;

% Define time (100 time steps for now)
t = -0.5:0.01:2; 

% Lift function to be used for comparison 
L = 2*sqrt((U*t)./c.*(1-(U.*t./c)));

w = 0:50;
S = 1./sqrt(1+2*pi*w);
S2 = sqrt((1./(1+(2*pi*w).^1.3)).^(1/1.3));

F = ((real(fer2).^2+imag(fer2).^2).^0.5);


figure (1)
plot(t,L);
ylim([-0.2,2])
grid on

figure(2)
plot(w,S2,'x',f,F,'o')
set(gca,'FontSize',16)
ylabel('Transform of Lift')
xlabel('F')
grid on
legend('Sears Function','BVI Data')


% %% Validation Method as Presented in Unsteady Blade Response:The BVI Model vs. the Gust Model (Grace, 2000)
% 
% % Analytical Lift Response of a Flat Plate to a passing vortex
% 
% % Establish normalized frequency k = (wc/2U) vector
% k = 0:1:50;
% 
% % Establish Parameters 
% gamma = 0.02; % Strength of Vortex
% c = 1; %chord length
% U = 1; % Freestream
% rho = 1.225; % Density of air
% h = 0.04; % Closest approach distance of vortex
% 
% H0 = besselh(0,1,k);
% H1 = besselh(1,1,k);
% S = 2./(pi*k.*(H0.^2-1i.*H1.^2));
% L = rho*gamma*(c/2)*exp((-k*h)).*(-1i.*S);
% 
% figure(100)
% plot(k,L)
% 
% 
% 




