
clear 
close all
load L_unsteady_9x
D1 = BEM2d_readbin('2NACA0001_long');
phi_dot_calc;
figure(25)
plot(L_unsteady);
hold on;
plot(D1.L);

test_low        = ones(1,100);
test_high       = ones(1,100);
test_low(50)    = 0.75;
test_high(50)   = 0.65;

figure(11)
hp1 = plot(test_low,'r');
hold on;
hp2 = plot(test_high,'k');
set(gca,'FontSize',15);
set(hp1,'LineWidth',2);
set(hp2,'LineWidth',2);
grid on;
ylim([0.6 1.1]);
legend('Low Signal','High Signal');

% high signal
time = 0.01:0.01:1;
foft = test_high;
ffer3;
fer_high = fer2;
f_high = f;

foft = test_low;
ffer3;
close figure 1
fer_low = fer2;
f_low = f;

figure(22)
hp1 = semilogy(f_high,abs(fer_high),'k');
hold on;
hp2 = semilogy(f_low,abs(fer_low),'r');
grid on;
set(gca,'FontSize',16);
xlabel('Frequency - f');
ylabel('FFT');
leg1 = legend('High Signal','Low Signal');