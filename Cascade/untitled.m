N = 10;
x = [-2:0.15:2];
y = x;
gam = 1;


xmat =[]
for i = 1:length(y)
    xmat = [xmat; x];
end

ymat = transpose(xmat)

u = gam/2/pi*( ymat)./(xmat.^2 + ymat.^2);
v = - gam/2/pi*xmat./(xmat.^2 + ymat.^2);

figure(88)
quiver(x,y,u,v)
