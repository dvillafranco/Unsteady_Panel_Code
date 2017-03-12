
Panel_BEM_Sears

clear L_unsteady

load Panel_Unsteady_noconvect2.mat


phi_dot_calc;
Lift_Panel_Interp;

figure(2000)
plot(-L_unsteady_fix);

figure(10)
plot(time,-L_unsteady_fix,'LineWidth',2);
legend('dt = 0.005','BEM','dt = 0.01');

clear foft
foft = -L_unsteady_fix;
