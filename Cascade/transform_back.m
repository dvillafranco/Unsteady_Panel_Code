

z2_shape = xnew+ 1i*ynew;
z1_shape = h/pi*atanh(z2_shape);
x_shape = real(z1_shape);
y_shape = imag(z1_shape);

z2_fixed = x_fixed+1i*y_fixed;
z1_fixed = h/pi*atanh(z2_fixed);
x_fx = real(z1_fixed);
y_fx = imag(z1_fixed);

y_fx = y_fx(1:end-2);
x_fx = x_fx(1:end-2);


figure(1)
plot(x_shape,y_shape)
hold on
plot(x_fx,y_fx)

z2_imposed = x_imp + 1i*y_imp;
z1_imposed = h/pi*atanh(z2_imposed);
x_imposed = real(z1_imposed);
y_imposed = imag(z1_imposed);
plot(x_imposed,y_imposed)
grid on
axis equal

title('Z1 space')




figure(2)

hold on
plot(x_fixed,y_fixed,'*')
plot(x_imp,y_imp,'o')
plot(xnew,ynew)
title('Z2 Space')
grid on
axis equal


