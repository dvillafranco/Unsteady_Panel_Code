% This script will comapre the fft plots for the NACA 0012 airfoil with the
% correction FFT from lysack? 
Panel_BEM_Sears;
k = 0.01:0.01:100;
 % Thickness Correction 
t_max = 0.12;
%S = (1./(1+(2*pi*k).^1.3)).^(1/1.3);
x_c = [0.00,0.000001];
tau_c = zeros(1,2);
for i= 1:2;
tau_c(i) = t_max/0.2*(0.2969*(x_c(i)^0.5)-0.1260*(x_c(i)^2)+0.2843*(x_c(i)^3)-...
    0.1015*(x_c(i)^4));
end

D_sum = sum(tau_c);

beta_atten = (2*pi*log(10)/10)*(1.123*D_sum+5.317*D_sum^2);

H_resp = (abs(S)).^2.*exp(-1*beta_atten*k/pi).*(1+0.8*t_max)^2;
H_mag = (real(H_resp).^2+imag(H_resp).^2).^0.5;

% figure(30)
% 
% hp1 = semilogy(k,abs(H_resp)*pi);
% xlabel('Reduced Frequency');
% ylabel('Magnitude of FFT')
% grid on;
% set(gca,'FontSize',16);

% 
% 
% 
% hp1 = semilogy(f_Panel*pi,abs(fer_Panel/vortstrength)./exp(-f*h*pi/1)*pi);
% 

