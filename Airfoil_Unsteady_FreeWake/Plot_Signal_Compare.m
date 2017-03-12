load L_unsteady_9x_airfoildensity.mat

phi_dot_calc;
time = 0.01:0.01:20;
D1 = BEM2d_readbin('2NACA0001_long');
figure(1)
hp1 = plot(D1.t,D1.L);
hold on;
hp2 = plot(time, L_unsteady);
set(hp1,'LineWidth',2);
set(hp2,'LineWidth',2);
xlabel('Time(s)')
ylabel('Lift coefficient')
grid on;
set(gca,'FontSize',16)
leg1 = legend('BEM Signal','Panel Code Signal');



