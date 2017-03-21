
dt = 0.01;

x_imp(1) = -9.0;
y_imp(1) = 0.02;


for gt = 1:4000;
    xp =  x_imp(gt);
    yp = y_imp(gt);
    psiinfield_vort_wake;
    x_imp(gt+1) = x_imp(gt) + upan*dt;
    y_imp(gt+1) = y_imp(gt) + vpan*dt;
end

x_fixed(1) = xnew(end) +0.45*dt;
y_fixed(1) = 0;

for gt = 1:4000;
    xp =  x_fixed(gt);
    yp = y_fixed(gt);
    psiinfield_vort_wake;
    x_fixed(gt+1) = x_fixed(gt) + upan*dt;
    y_fixed(gt+1) = y_fixed(gt) + vpan*dt;
end


% This will take the determined paths in the airfoil space and conduct the
% transform of the paths into the cascade plane.

% Establish complex path (x,) of imposed vortex & wake vorticies
vort_imp_z1 = x_imp + 1i*y_imp;
vort_wake_z1 = x_fixed + 1i*y_fixed;

% Imposed vortex path in the transformed (z2) plane
vort_imp_z2 = tanh(pi * vort_imp_z1/h);
% Wake vorticies path in the transformed (z2) plane
vort_wake_z2 = tanh(pi *vort_wake_z1/h);


figure(1000)
plot(xnew,ynew);
hold on;
plot(real(vort_imp_z2),imag(vort_imp_z2));
plot(real(vort_wake_z2),imag(vort_wake_z2));
plot(real(vort_imp_z1),imag(vort_imp_z1));
plot(real(vort_wake_z1),imag(vort_wake_z1));

