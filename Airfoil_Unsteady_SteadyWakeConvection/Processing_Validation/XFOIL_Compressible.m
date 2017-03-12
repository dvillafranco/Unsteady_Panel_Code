% Validation for NACA 0012 - Compressible (Single Airfoil Case)

% Define X-Foil Data for Alpha = 2.0 degrees
Mach_X = 0.05:0.05:0.6;
Cl_X = [0.0,0.2431, 0.2451, 0.2480, 0.2519, 0.2568, 0.2630, 0.2707,...
    0.2802, 0.2920, 0.3068, 0.3256];
Cl_BVI = [0.2413, 0.2423, 0.2440, 0.2464, 0.2497, 0.2538, 0.2589,...
    0.2652,0.2729, 0.2823,0.2939,0.3083];

error = abs(Cl_X-Cl_BVI)./Cl_X*100;
lbl_dwn = .1*max(Cl_X);
    
figure(1)
plot(Mach_X, Cl_X,'-o',Mach_X,Cl_BVI,'-*')
for i = 1:length(Mach_X)
    hold on
    %plot(Mach_X(i), Cl_X(i),'-o',Mach_X(i),Cl_BVI(i),'-*')
    %text(Mach_X(i)+0.008, Cl_BVI(i),num2str(error(i)));
    
    %text(Mach_X(i)+0.008, Cl_BVI(i),num2str(error(i)));
    xlabel('Mach Number')
    ylabel('Coefficient of Lift')
    legend('X Foil Data','BVI Data')
    grid on
    axis equal
    set(gca,'fontsize',14)
    
end

figure(2)
% Define Xfoil data for various angles of attack (Mach 0.2 for various
% alpha)
Alfa = 0:2:10;
Cl_Xalfa = [0, 0.2480, 0.4964, 0.7454, 0.9955, 1.2474];
Cl_BVIalfa = [1.3221e-15,0.2464, 0.4917, 0.7345, 0.9738, 1.2083];
figure(2)
plot(Alfa,Cl_Xalfa,'-o',Alfa,Cl_BVIalfa,'-*')
xlabel('Angle of Attack (Degrees)')
ylabel('Coefficient of Lift')
grid on
legend('X Foil Data','BVI Data')
set(gca,'fontsize',14)




