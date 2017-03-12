
h = 1;

zshape  = xnew + 1i*ynew;
z2shape = tanh(pi*zshape/h);
x2shape = real(z2shape);
y2shape = imag(z2shape);

zimp = x_imp + 1i*y_imp;
z2imp = tanh(pi*zimp/h);
x2imp = real(z2imp);
y2imp = imag(z2imp);

zfixed = x_fixed + 1i*y_fixed;
z2fixed = tanh(pi*zfixed/h);
x2fixed = real(z2fixed);
y2fixed = imag(z2fixed);


figure(1)
plot(xnew,ynew)
hold on
plot(x_imp,y_imp)
plot(x_fixed,y_fixed)
axis equal

figure(2)
plot(x2shape,y2shape)
hold on 
plot(x2imp,y2imp)
plot(x2fixed,y2fixed,'*')
axis equal
title('Z2 Space')