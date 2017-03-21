clear xnew ynew thetpan thetpan1 A B D C E R T11 T2 Gd1 Gd2 f4 f3 xmidmat
%h = 10;
% get airfoil points to build the cascade
%

% NACA 4 digit series for now 


t = .01;
% m = .02;
% p = .4;
m = .00;
p = .0;

%chi = -0*pi/180;

%vinf = 1;  % freestream velocity
%alp = 10*pi/180;   % freestream angle of attack

% 
% dx = .01;
% 
% xf = [0:dx:1*p];
% xa = [1*p+dx:dx:1];
% 
% x = [xf xa];

% this will work for symmetric airfoils 
% theven = [0 : .04 : pi];
% x = - cos(theven)/2 + .5;
% xf = 0;
% xa = x(2:end); 

%theta for break point  p = -cos(theven)/2 + .5

%dth = 0.04
xnums = 100;

thbreak = acos(-(p - .5)*2);

thfrac = thbreak/pi;
nbefore = round(xnums*thfrac)
nafter = xnums - nbefore;
if nbefore == 0 
    dthbef = 0
else
dthbef = thbreak/(nbefore);
end
thf = [0 : dthbef : thbreak];
dthaft = (pi - thbreak)/(nafter);
tha = [thbreak + dthaft : dthaft : pi];


%  xf = [0:.001:p];
%  xa = [p+.001:.001:1];

 xf = -cos(thf)/2 + .5;
 xa = -cos(tha)/2 + .5;
x = [xf xa];
numb = length(x);


yt = t/.2*(.2969*sqrt(x) - .126*x - .35160*x.^2 + .28430*x.^3 - .1015*x.^4);
if (m == 0) 
    ycf = 0*xf;
    yca = 0*xa;
    dycf = 0*xf;
    dyca = 0*xa;
else
    ycf = m/p^2 *(2*p*xf - xf.^2);
    yca = m/(1-p)^2*( (1-2*p) + 2*p*xa - xa.^2);


    dycf = m/p^2*(2*p-2*xf);
    dyca = m/(1-p)^2*(2*p - 2*xa);
end

yc = [ycf yca];
dyc = [dycf dyca];

thet = atan(dyc);


xu = x - yt.*sin(thet);
yu = yc + yt.*cos(thet);

xl = x + yt.*sin(thet);
yl = yc - yt.*cos(thet);



% Now put points from trailing edge around bottom and then around top back
% to trailing edge. 

xaf = [fliplr(xl) xu(2:end)] -.5 ;
yaf = [fliplr(yl) yu(2:end)];

%close trailing edge
yave = (yaf(end) - yaf(1))/2 + yaf(1);
yaf(end) = yave;
yaf(1) = yave;
xaf(end) = xaf(1);
% Add stagger if there is any
%chi = -00*(pi/180);

xaf2 = xaf*cos(chi) - yaf*sin(chi) ;
yaf2 = xaf*sin(chi) + yaf*cos(chi) ;


%this is a fudge here so that the asymp solver works - it needs airfoil
%going from 0 to 2 even after it is rotated (I think ) 
stretchfac = (max(xaf2)-min(xaf2))/1.0;
xaf2s = xaf2/stretchfac;
yaf2s = yaf2/stretchfac;

SJ_af = sqrt( ( xaf2(2:end) - xaf2(1:end-1) ).^2 + (yaf2(2:end)-yaf2(1:end-1)).^2);

figure(202)
plot(xaf2s,yaf2s)
%title('xaf2s, yaf2s','FontSize',16)
title('Cascade of Airfoils','FontSize',16)
%axis equal;
 hold on;
% set spacing (based on chord = 1) and repeat to create cascade of foils

%h = 1;
%stop
for j = 1:4
    yafs(j,:) = yaf2s +(j -1)*h;
    
    plot(xaf2s, (yafs(j,:)));
   
end
xnew = xaf2;
ynew = yaf2;

% set trailing edge point 
slope2 = (ynew(3) - ynew(2))/(xnew(3) - xnew(2));
slope3 = (ynew(end-1) - ynew(end-2))/(xnew(end-1) - xnew(end-2));
xnew(1) = (-ynew(2) + ynew(end-1) + xnew(2)*slope2 - xnew(end-1)*slope3)/ ...
   (slope2 - slope3);
ynew(1) = slope2 * (xnew(1) - xnew(2)) + ynew(2);

xnew(end) = xnew(1);
ynew(end) = ynew(1);
set(gca,'FontSize',16)
% Form the complex variable and do the transformation

z1 = xnew + 1i*ynew;

z2 = tanh(pi * z1/h);
xnew = real(z2);
ynew = imag(z2);
figure(1)
%plot(xnew,ynew,'+');
plot(xnew,ynew);
axis equal
title('After transformation','FontSize',16)
set(gca,'FontSize',16)
  

% Determine the panel angles used in Keuthe and Chow 

for j = 1: length(xnew) -1;
    
        thetpan1(j) = atan2(ynew(j+1) - ynew(j), xnew(j+1) - xnew(j));
end;
thetpan = (thetpan1);
figure(3)
plot(thetpan*180/pi,'+')

% xnew and ynew and thetpan are the end points and panel angles.
% from this must get panel midpoints and lengths

xmid = xnew(1:end-1) + (xnew(2:end) - xnew(1:end-1))/2;
ymid = ynew(1:end-1) +  (ynew(2:end) - ynew(1:end-1))/2;


SJ = sqrt( ( xnew(2:end) - xnew(1:end-1) ).^2 + (ynew(2:end)-ynew(1:end-1)).^2);



