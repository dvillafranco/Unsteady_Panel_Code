%%%%%%%%%%%%%%%%%%%%%%%%%%Ventres_k2_h.m%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Date: 03/17/16
% Version: 1.00
%
% Author: Dorien O. Villafranco
% Unsteady Fluid Mechanics and Acoustics Group
% Department of Mechanical Engineering 
% College of Engineering
% Boston University
% 
% Program Description: 
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
close all


%% Variation of K1 and h
% Load the Ventres File of interest
filename    = 'Ventres_k130_h10.txt';               % File name for data
A           = importdata(filename);                 % Load data
h_num       = 10;                                   % # of spacings (h)
k1_num      = 30;                                   % # of k1's 
cl_ven      = zeros(k1_num,h_num);                  % Initialize cl_ven
k_ven       = linspace(1,k1_num,k1_num);            % Initialize k_ven

cl_ven = reshape(A(:,2),size(cl_ven));
figure(1)
for i = 1:h_num
plot(k_ven,cl_ven(:,i),'LineWidth',2)
hold on;
end
legend('h=1','h=2','h=3','h=4','h=5','h=6','h=7','h=8','h=9','h=10');
set(gca,'FontSize',16);
grid on;
xlabel('K1')
ylabel('Ventres Lift Coefficient');


% Set up sigma formulation
stag = 0;
k1 = 1;
k2 = 1:30;
h = 1:10;
for j = 1:length(h)
    h_chalf = 2*h(j);
    for i = 1:length(k2)
    sigma(i,j) = k1*h_chalf*sin(stag) + k2(i)*h_chalf*cos(stag);
    end
end
figure(3)
plot(k2,sigma,'LineWidth',2)
hold on;
plot(k2,100*ones(size(k2)),'k','LineWidth',2)
grid on
set(gca,'FontSize',16)
xlabel('k2')
ylabel('sigma')
legend('h=1','h=2','h=3','h=4','h=5','h=6','h=7','h=8','h=9','h=10','Location','NorthWest');
for i = 1:10
[xout,yout] = intersections(k2,sigma(:,i),k2,100*ones(size(k2)),1);
plot(xout,yout,'r.','markersize',25)
end




%% Variation of K2 and h
% Load the Ventres File of interest
% Load the Ventres File of interest
filename    = 'Ventres_k230_h10_2.txt';               % File name for data
A           = importdata(filename);                 % Load data
h_num       = 10;                                   % # of spacings (h)
k2_num      = 290;                                  % # of k2's 
cl_ven      = zeros(k2_num,h_num);                  % Initialize cl_ven
k_ven       = linspace(1,k1_num,290);               % Initialize k_ven

cl_ven = reshape(A(:,2),size(cl_ven));
% figure(2)
% for i = 1:h_num
% plot(k_ven,cl_ven(:,i),'LineWidth',2)
% hold on;
% end
% legend('h=1','h=2','h=3','h=4','h=5','h=6','h=7','h=8','h=9','h=10');
% set(gca,'FontSize',16);
% grid on;
% xlabel('K2')
% ylabel('Ventres Lift Coefficient');

figure(5)
plot(k_ven,cl_ven(:,2),'LineWidth',2)
hold on;
plot(k_ven,cl_ven(:,5),'LineWidth',2)
y1=get(gca,'ylim');
plot([10,10],y1,'k','LineWidth',2);
plot([25,25],y1,'k','LineWidth',2);
xlabel('K2');
ylabel('Ventress Lift')
legend('h=2','h=5')
grid on
set(gca,'FontSize',16)

