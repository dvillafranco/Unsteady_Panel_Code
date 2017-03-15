% Unsteady Coefficient of Pressure

% This script calculates the unsteady coefficient of pressure for the panel
% code to calculate unsteady lift for the blade vortex interaction problem 
dt = 0.01;
% Phi on the surface of the body is stored in time along with the
% velocities associated with the vortex singularities defining the airfoil.
         z2 = xpmat + 1i*ypmat;
         z1 = atanh(z2)*h/pi;
         %z1 = xpinaf + 1i*ypinaf;
         dz2dz1 = pi/h*(1 - (tanh(pi*z1/h)).^2 );
         
         compv2 = upan_body_time - 1i*vpan_body_time;
%          compv1 = compv2.*dz2dz1;
         compv1 = bsxfun(@times,compv2,dz2dz1);
         
         u1 = real(compv1);
         v1 = -imag(compv1);

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
% v_test_total = (upan_body_time.^2 + vpan_body_time.^2);
v_test_total = u1.^2 + v1.^2;
    

%Cp_unsteady = 1 - 2*phi_dot - (upan_body_time.^2 + vpan_body_time.^2);
Cp_unsteady = 1 - 2*phi_dot(:,1:end) - v_test_total;
for ipt = 1:size(Cp_unsteady,2)
    L_unsteady(ipt) = sum(Cp_unsteady(:,ipt)'.*SJ_af.*cos(thetpan))/vinf^2;
end