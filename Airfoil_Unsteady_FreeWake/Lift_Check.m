% Plot check 
clear all
close all


load Unsteady_NACA0001.mat
time = 0.0125:0.0125:10;

lift_gam = zeros(800,1);
for i = 1:length(lift_coef)
    lift_gam(i) = 2*kelvin(1:end-1)*gamans_MAT(1:end-1,i);
end

figure(1)
    plot(time,-lift_gam,'x','LineWidth',2)
    hold on;
    plot(time,lift_coef,'LineWidth',2)
    xlabel('Time (s)')
    ylabel('Coefficient of Lift')
    set(gca,'FontSize',15)
    
    grid on;
    
D = BEM2d_readbin('2example');
plot(D.t,D.L,'--','LineWidth',2);


legend('Cl_{Gamma}','Cl_{Pressure}','Cl_{BEM}','Location','southwest')