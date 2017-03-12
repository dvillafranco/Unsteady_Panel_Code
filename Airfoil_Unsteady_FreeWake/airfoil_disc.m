
% get airfoil points to build the cascade
%

% NACA 4 digit series for now 

t = .01;
m = .00;
p = .0;
% m = .11;
% p = .45;
chi = 0;

%vinf = 1;  % freestream velocity
%alp = 0*pi/180;   % freestream angle of attack


%chi = -20*pi/180;
%chi = -6*pi/180;
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
xnums = 80;


thbreak = acos(-(p - .5)*2);

if thbreak == 0
    thbreak=1;
end


thfrac = thbreak/pi;
nbefore = round(xnums*thfrac);
nafter = xnums - nbefore;
dthbef = thbreak/(nbefore);
thf = [0 : dthbef : thbreak];
dthaft = (pi - thbreak)/(nafter);
tha = [thbreak + dthaft : dthaft : pi];
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

plot(xu,yu, xl,yl);
axis equal;

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


% put in the PG transform here ... (then the chi)
%yaf = yaf/beta;
%xaf = xaf/beta;

xaf2 = xaf*cos(chi) - yaf*sin(chi) ;
yaf2 = xaf*sin(chi) + yaf*cos(chi) ;

figure(1)
plot(xaf2,yaf2)
%axis equal;
 hold on;
% set spacing (based on chord = 1) and repeat to create cascade of foils


% h = 1;
% for j = 1:4
%     yafs(j,:) = yaf2 +(j -1)*h;
%     
%     plot(xaf2, (yafs(j,:)));
%    
% end
% 
% % Form the complex variable and do the transformation
% 
% z1 = xaf2 + 1i*yaf2;
% 
% z2 = tanh(pi * z1/h);
% xnew = real(z2);
% ynew = imag(z2);
% figure(2)
% plot(xnew,ynew,'+');
% axis equal
  

xnew = xaf2;
ynew = yaf2;

% set trailing edge point 
% slope2 = (ynew(3) - ynew(2))/(xnew(3) - xnew(2));
% slope3 = (ynew(end-1) - ynew(end-2))/(xnew(end-1) - xnew(end-2));
% xnew(1) = (-ynew(2) + ynew(end-1) + xnew(2)*slope2 - xnew(end-1)*slope3)/ ...
%    (slope2 - slope3);
% ynew(1) = slope2 * (xnew(1) - xnew(2)) + ynew(2);
% 
% xnew(end) = xnew(1);
% ynew(end) = ynew(1);
% % 
% % figure(20)
% plot(xnew,ynew,'--x');

spacing = xnew(52) - xnew(51); 

xnew  = [0.5:-0.04:-0.18 , xnew(52:110), -0.18:0.04:0.5];
ynew  = [0: -0.000294: -0.004998 , ynew(52:110), 0.004998:-0.000294:0,];


xtemp = xnew;
ytemp = ynew;



% xtemp = [ xnew(1) xnew(14) xnew(21) xnew(22:end-21) xnew(end-20) xnew(end-13) xnew(end) ] ;
% ytemp = [ ynew(1) ynew(14) ynew(21) ynew(22:end-21) ynew(end-20) ynew(end-13) ynew(end) ] ;
%this is a fudge here so that the asymp solver works - it needs airfoil
%going from 0 to 2 even after it is rotated (I think ) 
stretchfac = (max(xtemp)-min(xtemp))/1.0;
xnew = xtemp/stretchfac;
ynew = ytemp/stretchfac;


% ynew(1) = 0;
% ynew(end) = 0;
% xnew(1) = .505;
% xnew(end) = .505;
% Determine the panel angles used in Keuthe and Chow 

for j = 1: length(xtemp) -1;
    
        thetpan1(j) = atan2(ytemp(j+1) - ytemp(j), xtemp(j+1) - xtemp(j));
end;
thetpan = unwrap(thetpan1);
figure(3)
plot(thetpan*180/pi,'+')

% xnew and ynew and thetpan are the end points and panel angles.
% from this must get panel midpoints and lengths

xmid = xtemp(1:end-1)+ (xtemp(2:end) - xtemp(1:end-1))/2;
ymid = ytemp(1:end-1)+(ytemp(2:end) - ytemp(1:end-1))/2;

SJ = sqrt( ( xtemp(2:end) - xtemp(1:end-1) ).^2 + (ytemp(2:end)-ytemp(1:end-1)).^2);

clear xnew ynew 
xnew = xtemp;
ynew = ytemp;

fileID = fopen('NACA0001_mixed.dat','w');
for i = 1:length(xnew)
    fprintf(fileID,'%f %f\n',xnew(i),ynew(i));
end
fclose(fileID);

figure(100)
plot(xnew,ynew,'x')

