% Dorien Villafranco
% Boston University 
% Department of Mechanical Engineering

% This script will plot points of airfoil xpmat, ypmat along with path of
% vortex travelled.

% Another subplot will be introduced of lift values along the airfoil

xplot = linspace(-0.5,0.5,1805);

yimp_plot = 0.1*ones(size(x_imp));
subplot(2,1,1)
plot(xnew,ynew,x_imp,yimp_plot,'LineWidth',1.5)
ylim([-0.1,0.1])
xlim([-3,2])
xlabel('x')
ylabel('y')
title('Path of Vortex over Airfoil')
%axes equal

set(gca,'FontSize',16)

% subplot(2,1,2)
% plot(time,foft,'LineWidth',1.5)
% xlim([0,2.5])
% xlabel('Time')
% ylabel('Lift')
% set(gca,'FontSize',16)