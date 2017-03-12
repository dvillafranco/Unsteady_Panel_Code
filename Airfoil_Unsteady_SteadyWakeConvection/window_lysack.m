% This will take the response of lift_coef from main script and do the 
% following:
% 1. Create the window function as shown in Lysack, 2013
% 2. Multiply the window function by the impulse function to create the 
% tapered function
% 3. Plot the impulse response, window function, and tapered repsponse 
% 4. Take the fft of the tapered function and plot results 
% Plots should be similar to figures 5 and 6 in Lysack, 2013. 

time_disc = linspace(0,10.01,1001);
gamma = gam_imp;
t = 4.5:0.01:5.5;
rho = 1;
h = 0.02;
k = 0.01:0.1:100.01;
w_bef = zeros(1,350);
w_aft = zeros(1,550);
U = 1; 
c = 1;
w_fxn = [w_bef, (cos((pi*U*t)/c)).^2,w_aft];
L = lift_coef*0.5;
L = L/(gamma) ;


% Sears Forumlation
H0 = besselh(0,2,k);
H1 = besselh(1,2,k);
S = 2./(pi*k.*(H0-(1i*(H1))));
%S = ((2*i)./(pi*k)).*(exp(-i*k)./(H1+i*H0));
L_2 = rho*(c/2)*(-1i*S).*exp(-k*h);
L_mag = abs(L_2);



taper_fxn = L.*w_fxn;
%taper_fxn = L;

figure(1)
plot(time_disc,taper_fxn,'-','LineWidth',2)
hold on 
plot(time_disc,L,'--','LineWidth',2)
plot(time_disc,w_fxn,'LineWidth',2)
legend('Taperd Response','Impulse Response','Window Function')
grid on

figure(3)
plot(time_disc,taper_fxn)

L_fft = fft(taper_fxn);
L_fft_mag = abs(L_fft)/(2*pi);


figure(2)
f = (0:length(L_fft_mag)-1)*1001/length(L_fft_mag);
semilogy(f*pi,L_fft_mag,k,L_mag)
xlim([0.01 100])
%ylim([10^-2, 1])
grid on