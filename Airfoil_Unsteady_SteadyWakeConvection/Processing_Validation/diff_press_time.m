% this script will plot the pressure surface data for a NACA 0012 airfoil
% in the unsteady conditions. 

load Run_Data/NACA0012_Pressure_Data.mat

%cp_time_diff = cp_time(1:120,:) - cp_time(121:end,:);
gamma = 0.01;

for q = 1:240
    foft = cp_time(q,:);
    ffer
    cp_trans(q,:) = ((real(fer2)/gamma).^2 + (imag(fer2)/gamma).^2).^0.5;
end

cp_freq_diff = -1* cp_trans(1:120,:) + cp_trans(end:-1:121,:);


figure(2222)
plot(cp_freq_diff(:,500));

stop




h = 0.12;

k = linspace(0,60,1100);


F = ((real(fer2/gamma).^2+imag(fer2/gamma).^2).^0.5).*exp(-k*h);

semilogy(k,F)