clear;
close all


% in order to run a compressible flow case 
% set the right NACA or whatever in casc_justaf.mtemp1

sos = 340;
mach = 0.0;
beta = sqrt(1 - mach^2);

if (mach == 0.0) 
       vinf = 1.0;
else
    vinf = mach*sos;
    %vinf = 1.0;
end

alp = 0*pi/180;   % freestream angle of attack

%alp = -atan((-0.02/(2*pi*4.5)))

%chi = -20*pi/180;
chi = -0*pi/180;
%chi = -13.6*pi/180


% get the geometry and readjust the panel information
philocnew_store = zeros(120,800);
casc_justaf;
%airfoil_disc;
figure(41)
plot(xnew,ynew,'x')
hold on;

%ynew = ynew*beta;
% xnew = xnew*beta;
% 
plot(xnew,ynew,'r')
% Detervmine the panel angles used in Keuthe and Chow 

for j = 1: length(xnew) -1;
    
        thetpan1(j) = atan2(ynew(j+1) - ynew(j), xnew(j+1) - xnew(j));
end;
thetpan = unwrap(thetpan1);
figure(3)
plot(thetpan*180/pi,'+')

% xnew and ynew and thetpan are the end points and panel angles.
% from this must get panel midpoints and lengths

xmid = xnew(1:end-1)+ (xnew(2:end) - xnew(1:end-1))/2;
ymid = ynew(1:end-1)+(ynew(2:end) - ynew(1:end-1))/2;

SJ = sqrt( ( xnew(2:end) - xnew(1:end-1) ).^2 + (ynew(2:end)-ynew(1:end-1)).^2);



% get the values of gamma on the panels

find_gams_justaf_only;

psionsurface_original;

%find_gams_justaf_withvortex;


%Wagner_Test;


% find phi and psi on the surface

%psionsurface;

% xpmat = xmidmat(1:1:end,1) -.000001*sini(1:1:end,1);
% ypmat = ymidmat(1:1:end,1)  + .000001*cosi(1:1:end,1);

stop
%compute lift in actual plane
cp_comp = cp/beta;
yaf_comp = yaf*beta;
xaf_comp = xaf;


xaf2_comp = xaf_comp*cos(chi) - yaf_comp*sin(chi) ;
yaf2_comp = xaf_comp*sin(chi) + yaf_comp*cos(chi) ;

for j = 1: length(xnew) -1;
        thetpan1_comp(j) = atan2(yaf2_comp(j+1) - yaf2_comp(j), xaf2_comp(j+1) - xaf2_comp(j));
end;
thetpan_comp = unwrap(thetpan1_comp);

SJ_comp = sqrt( ( xaf2_comp(2:end) - xaf2_comp(1:end-1) ).^2 + (yaf2_comp(2:end)-yaf2_comp(1:end-1)).^2);
lift_coef_com = cp_comp'.*SJ_comp.*cos(thetpan_comp)/vinf^2;
lift_coef_comp = sum(cp_comp'.*SJ_comp.*cos(thetpan_comp))/vinf^2;

stop
 % now xpmat, ypmat are the locations on the surface of the airfoil in the
 % PG plane (incompressible plane) Find incompressible F^i from Meyers
 % appendix
        
        
        
    
        psibeta = psiloc/beta^2;
        phibeta = philocnew/beta;
        
         phi1i = (philocnew - vinf*xpmat)/vinf;
         psi1i = (psiloc - vinf*ypmat)/vinf;
         
         phi1 = phi1i/beta;
         psi1 = psi1i/beta^2;
         
%         figure(35)
%         hold on;
        
        % the leading edge point is a node so the collocation point - where
        % the values are computed will never be at the leading edge exactly
        % unless the nodes are shifted somehow
        [val leloc] = min(xpmat);
        
        plot((xpmat+.5)*2, (phi1-phi1(leloc))*2,'b', (xpmat+.5)*2, (psi1-psi1(leloc))*2,'r');  % need multiply by 2 to base on chord of 2
       
      %Must unrotate, take out of PG plane and rotate again
      xunrot = xpmat*cos(-chi) - ypmat*sin(-chi) ;
      yunrot = xpmat*sin(-chi) + ypmat*cos(-chi) ;

      yunrot = yunrot/beta;
      xfin = xunrot*cos(chi) - yunrot*sin(chi) ;
        yfin = xunrot*sin(chi) + yunrot*cos(chi) ;
% 
%          compmat = [ (xfin+.5)*2 (yfin)*2 (phi1-phi1(leloc))*2  (psi1-psi1(leloc))*2 ];
%       %  dlmwrite('phi1_psi1.txt',compmat,'delimiter','\t','precision',8)
%       
%         compmattop = [ (xfin(leloc:end)+.5)*2, yfin(leloc:end)*2  (phi1(leloc:end)-phi1(leloc))*2 ]
%        dlmwrite('phitop_naca0006.txt',compmattop,'delimiter','\t','precision',8)
%       
%       compmatbot = [ flipud(xfin(1:leloc)+.5)*2, flipud(yfin(1:leloc))*2  flipud(phi1(1:leloc)-phi1(leloc))*2 ]
%        dlmwrite('phibot_naca0006.txt',compmatbot,'delimiter','\t','precision',8)