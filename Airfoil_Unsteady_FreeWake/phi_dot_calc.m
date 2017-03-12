% Unsteady Coefficient of Pressure

% This script calculates the unsteady coefficient of pressure for the panel
% code to calculate unsteady lift for the blade vortex interaction problem 

% Phi on the surface of the body is stored in time along with the
% velocities associated with the vortex singularities defining the airfoil.


% phi_dot = 0*phiII_overS_time;
% phi_dot(:,2:end-1) = phiII_overS_time(:,3:end) - ...
%     phiII_overS_time(:,1:end-2);
% phi_dot(:,1) = -phiII_overS_time(:,3) + ...
%     4*phiII_overS_time(:,2) - 3*phiII_overS_time(:,1);
% phi_dot(:,end) = phiII_overS_time(:,end-2) - ...
%     4*phiII_overS_time(:,end-1) + 3*phiII_overS_time(:,end);
% phi_dot = phi_dot/(2*dt);
% 
% Cp_unsteady = 1 - 2*phi_dot - (upan_body_time.^2 + vpan_body_time.^2);
% for ipt = 1:size(upan_body_time,2)
%     L_unsteady(ipt) = sum(Cp_unsteady(:,ipt)'.*SJ.*cos(thetpan))/vinf^2;
% end

% for i = 1:length(D)
%       D(i).phidot = 0*D(i).Pb;
%    D(i).phidot(:,2:end-1) = ( D(i).Pb(:,3:end) -   D(i).Pb(:,1:end-2) );
%       D(i).phidot(:,   1   ) = (-D(i).Pb(:,  3  ) + 4*D(i).Pb(:,   2   ) - 3*D(i).Pb(:, 1 ) );
%       D(i).phidot(:,  end  ) = ( D(i).Pb(:,end-2) - 4*D(i).Pb(:, end-1 ) + 3*D(i).Pb(:,end) );
%       D(i).phidot            = D(i).phidot / (2*dt);
% end


% phi_dot = 0*ansphipan_time;
% phi_dot(:,2:end-1) = ansphipan_time(:,3:end) - ...
%     ansphipan_time(:,1:end-2);
% phi_dot(:,1) = -ansphipan_time(:,3) + ...
%     4*ansphipan_time(:,2) - 3*ansphipan_time(:,1);
% phi_dot(:,end) = ansphipan_time(:,end-2) - ...
%     4*ansphipan_time(:,end-1) + 3*ansphipan_time(:,end);
% phi_dot = phi_dot/(2*dt);
% ansphipan_diff = (ansphipan_time(:,2:end) - ansphipan_time(:,1:end-1));
% 
% v_total = bsxfun(@rdivide,ansphipan_diff,SJ');
%     
% 
% %Cp_unsteady = 1 - 2*phi_dot - (upan_body_time.^2 + vpan_body_time.^2);
% Cp_unsteady = 1 - 2*phi_dot(:,2:end) - v_total.^2;
% for ipt = 1:size(Cp_unsteady,2)
%     L_unsteady(ipt) = sum(-Cp_unsteady(:,ipt)'.*SJ.*cos(thetpan))/vinf^2;
% end


phi_dot = 0*philocpannew_time;
phi_dot(:,2:end-1) = philocpannew_time(:,3:end) - ...
    philocpannew_time(:,1:end-2);
phi_dot(:,1) = -philocpannew_time(:,3) + ...
    4*philocpannew_time(:,2) - 3*philocpannew_time(:,1);
phi_dot(:,end) = philocpannew_time(:,end-2) - ...
    4*philocpannew_time(:,end-1) + 3*philocpannew_time(:,end);
phi_dot = phi_dot/(2*dt);
ansphipan_diff = (philocpannew_time(:,2:end) - philocpannew_time(:,1:end-1));

v_total = bsxfun(@rdivide,ansphipan_diff,SJ');
v_test_total = (upan_body_time.^2 + vpan_body_time.^2);
    

%Cp_unsteady = 1 - 2*phi_dot - (upan_body_time.^2 + vpan_body_time.^2);
Cp_unsteady = 1 - 2*phi_dot(:,1:end) - v_test_total;
for ipt = 1:size(Cp_unsteady,2)
    L_unsteady(ipt) = sum(Cp_unsteady(:,ipt)'.*SJ.*cos(thetpan))/vinf^2;
end

% % Load the BEM data to compare
% D1 = BEM2d_readbin('2NACA0001_shift');
% %time = 0.01:0.01:20;
% time = 0.0125:0.0125:20;
% figure(1)
% hp1 = plot(D1.t,D1.L);
% hold on;
% hp2 = plot(time,L_unsteady);
% set(gca,'FontSize',16)
% grid on;
% set(hp1,'LineWidth',2)
% set(hp2,'LineWidth',2)
% leg1 = legend('BEM Signal','Panel Code Signal');
% xlabel('Time (s)')
% ylabel('Lift')