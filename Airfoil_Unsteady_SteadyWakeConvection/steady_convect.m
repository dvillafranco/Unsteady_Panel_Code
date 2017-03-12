
dt = 0.005;

x_imp(1) = -9.0;
y_imp(1) = 0.02;


for gt = 1:4000;
    xp =  x_imp(gt);
    yp = y_imp(gt);
    psiinfield_vort_wake;
    x_imp(gt+1) = x_imp(gt) + upan*dt;
    y_imp(gt+1) = y_imp(gt) + vpan*dt;
end

xnew(end) = 0.5;
x_fixed(1) = xnew(end) +0.504*dt;
y_fixed(1) = 0;

for gt = 1:4000;
    xp =  x_fixed(gt);
    yp = y_fixed(gt);
    psiinfield_vort_wake;
    x_fixed(gt+1) = x_fixed(gt) + upan*dt;
    y_fixed(gt+1) = y_fixed(gt) + vpan*dt;
end
