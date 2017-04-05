% Clear workspace 
clear;
close all


%% Set the free-stream parameters

sos = 340;                      % Speed of Sound
mach = 0;                       % Mach Number
beta = sqrt(1 - mach^2);        % Beta - Compressibility Factor

if (mach == 0.0) 
       vinf = 1.0;
else
    %vinf = mach*sos;
    vinf = 1.0;
end

alp =0*pi/180;                  % Freestream angle of attack

h = 0.5;                         % Spacing in cascade
chi = -0*pi/180;                % Stagger angle in cascade
%chi = -13.6*pi/180

%% Define Airfoil Points
% Set Airfoil Properties in this code. 
casc_justaf;

% Airfoil Panel Angles 
% Determine the panel angles used in Keuthe and Chow 
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

SJ = sqrt( ( xnew(2:end) - xnew(1:end-1) ).^2 + (ynew(2:end)-...
    ynew(1:end-1)).^2);

%% Steady Solve
% Solve for the gammas on the airfoil body - steady case
find_gams_justaf_only;

% Obtain potential on the surface and velocities on the surface
psionsurface_original;

%% Steady Convect
% Determine the path of the imposed vortex and wake vorticies in the
% airfoil plane and transform in to the cascade plane using conformal
% transormation
steady_convect_casc;



%% Cascade Points
% Set the same airfoil NACA #s here and following code will generate
% cascade of said NACA airfoils with spacing h defined above
casc_transform;

%% Cascade Steady Solve
% get the values of gamma on the panels for the cascade

find_gams_casc_only;

% find phi and psi on the surface for the cascade

find3_psionsurface;

gam_dim = zeros(1,4001);
find_gams_casc_vortex_shed;
stop
%compute lift in actual plane
cp_comp = cpcasc/beta;
yaf_comp = yaf/beta;



xaf2_comp = xaf*cos(chi) - yaf*sin(chi) ;
yaf2_comp = xaf*sin(chi) + yaf*cos(chi) ;

for j = 1: length(xnew) -1;
    
        thetpan1_comp(j) = atan2(yaf2_comp(j+1) - yaf2_comp(j), xaf2_comp(j+1) - xaf2_comp(j));
end;
thetpan_comp = unwrap(thetpan1_comp);

SJ_comp = sqrt( ( xaf2_comp(2:end) - xaf2_comp(1:end-1) ).^2 + (yaf2_comp(2:end)-yaf2_comp(1:end-1)).^2);

lift_coef_comp = sum(cp_comp'.*SJ_comp.*cos(thetpan_comp))/vinf^2;


stop
 % now xpmat, ypmat are the locations on the surface of the airfoil in the
 % PG plane (incompressible plane) Find incompressible F^i from Meyers
 % appendix
        
        
        
    
%         psibeta = psiloc/beta^2;
%         phibeta = philocnew/beta;
        
         phi1i = (philocnew - vinf*xpmat)/vinf;
         psi1i = (psiloc - vinf*ypmat)/vinf;
         
         phi1 = phi1i/beta;
         psi1 = psi1i/beta^2;
         
        figure(35)
        hold on;
        
        % the leading edge point is a node so the collocation point - where
        % the values are computed will never be at the leading edge exactly
        % unless the nodes are shifted somehow
        [val leloc] = min(xpmat);
        
        z2 = xpmat + 1i*ypmat;
        z1 = atanh(z2)*h/pi;
        
    %    plot((xpmat+.5)*2, (phi1-phi1(leloc))*2,'b', (xpmat+.5)*2, (psi1-psi1(leloc))*2,'r');  % need multiply by 2 to base on chord of 2
    
    plot((real(z1)+.5)*2, (phi1-phi1(leloc))*2,'b', (real(z1)+.5)*2, (psi1-psi1(leloc))*2,'r');  % need multiply by 2 to base on chord of 2
    
       
      %Must unrotate, take out of PG plane and rotate again
      xunrot = real(z1)*cos(-chi) - imag(z1)*sin(-chi) ;
      yunrot = real(z1)*sin(-chi) + imag(z1)*cos(-chi) ;

      yunrot = yunrot/beta;
      xfin = xunrot*cos(chi) - yunrot*sin(chi) ;
        yfin = xunrot*sin(chi) + yunrot*cos(chi) ;

         compmat = [ (xfin+.5)*2 (yfin)*2 (phi1-phi1(leloc))*2  (psi1-psi1(leloc))*2 ];
      %  dlmwrite('phi1_psi1.txt',compmat,'delimiter','\t','precision',8)
      
        compmattopcasc = [ (xfin(leloc:end)+.5)*2, yfin(leloc:end)*2  (phi1(leloc:end)-phi1(leloc))*2 ]
       dlmwrite('phitop_casc_naca0006stag30.txt',compmattopcasc,'delimiter','\t','precision',8)
      
      compmatbotcasc = [ flipud(xfin(1:leloc)+.5)*2, flipud(yfin(1:leloc))*2  flipud(phi1(1:leloc)-phi1(leloc))*2 ]
       dlmwrite('phibot_casc_naca0006stag30.txt',compmatbotcasc,'delimiter','\t','precision',8)