
%%%%%%% Incompressible NACA 0012


Alfa = 0:2:10;
Cl_X = [0, 0.2416, 0.4829, 0.7235, 0.9634, 1.2020];
Cl_BVI = [0, 0.2405, 0.4808,0.7168, 0.9523, 1.1816];
error012 = abs(Cl_BVI-Cl_X)./Cl_X*100;

figure(1)
plot(Alfa,Cl_X,'-*',Alfa,Cl_BVI,'-o')
xlabel('Angle of Attack (degrees)')
ylabel('Coefficient of Lift')
grid on
legend('X Foil Data','BVI Data')
set(gca,'fontsize',14)


%%%%%% Incompressible NACA 0001 - XFOIL

figure(2)
Cl01 = [0, 0.2229, 0.4453, .6664, 0.8858, 1.1028];
ClBVI01 = [0, 0.2204, 0.4397, 0.6569, 0.8709, 1.0806];

error1 = abs(Cl01 - ClBVI01) ./ Cl01 * 100;

plot(Alfa, Cl01,'-*', Alfa, ClBVI01, '-o')
xlabel('Angle of Attack (degrees)')
ylabel('Coefficient of Lift')
grid on
legend('X Foil Data','BVI Data')
set(gca,'fontsize',14)


%%%%%% Incompressible NACA 2412 - XFOIL

figure(3)
Cl2412 = [0.2554, 0.4968,0.7376,0.9775, 1.2162, 1.4534];
    
ClBVI2412 = [0.2624, 0.5031, 0.7420, 0.9780, 1.2098, 1.4365];

error2412 = abs(Cl01 - ClBVI01) ./ Cl01 * 100;

plot(Alfa, Cl2412,'-*', Alfa, ClBVI2412, '-o')
xlabel('Angle of Attack (degrees)')
ylabel('Coefficient of Lift')
grid on
legend('X Foil Data','BVI Data')
set(gca,'fontsize',14)
