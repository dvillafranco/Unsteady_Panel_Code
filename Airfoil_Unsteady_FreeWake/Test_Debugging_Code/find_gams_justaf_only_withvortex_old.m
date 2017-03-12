% find_gams_justaf_only_withvortex


% set up matrix to find unknown gamma values 
% use casc_tranform.m to get the endpoints, midpoints, and panel angles and
%   h which is the cascade gap

% Define Parameters for the Vortex
gam_imp = 0.01;
x_imp = -2.5;
y_imp = 0.04;

num = length(xnew);

for i = 1:num-1
        A(i,:) = -(xmid(i) - xnew(1:end-1)).*cos(thetpan) - (ymid(i) - ynew(1:end-1)).*sin(thetpan);
        B(i,:) = (xmid(i) - xnew(1:end-1)).^2 + (ymid(i) - ynew(1:end-1)).^2;
        D(i,:) = cos(thetpan(i) - thetpan);
        C(i,:) = -sin(thetpan(i) - thetpan);
        E(i,:) = (xmid(i) - xnew(1:end-1)).*sin(thetpan) - (ymid(i) - ynew(1:end-1)).*cos(thetpan);
        
        for j =1:length(E)
            if abs(E(i,j)) < 1e-12
                E(i,j) = -1e-12;
              %  G(i,j) = 0;
            end
        end
        
        R(i,:) = (xmid(i) - xnew(1:end-1))*cos(thetpan(i)) + (ymid(i) - ynew(1:end-1)).*sin(thetpan(i));
        T11(i,:) = -(xmid(i) - xnew(1:end-1))/2.*cos(thetpan(i) - 2*thetpan) ...
            + (ymid(i) - ynew(1:end-1))/2.*sin(thetpan(i)-2*thetpan) ;
        T2(i,:) = - (  (xmid(i) - xnew(1:end-1)).*sin(thetpan)-(ymid(i) - ynew(1:end-1)).*cos(thetpan) ).*cos(thetpan(i)-thetpan) ...
            +   (   (xmid(i) - xnew(1:end-1)).*cos(thetpan)-(ymid(i) - ynew(1:end-1)).*sin(thetpan)).*sin(thetpan(i)-thetpan);
end;


% create SJ matrix (repeat in each column)

SJmat = [];
for i = 1:num-1;
        SJmat = [SJmat;  SJ];
end;

t1 = 2*A.*SJmat+SJmat.^2 + B;

I1ij = -1/2*D.*(log(t1./B)) + (atan2(E.*SJmat,(B+SJmat.*A))).*C; 

t6 = (2*A.^2.*D + A.*R - B.*D)./E;
I2ij = 1./SJmat.*(-D.*SJmat + (log(t1./B)).*T11 - (atan2(E.*SJmat,(B+SJmat.*A))).*(A.*C - D.*E));

% get xmid , ymid and normal terms in matrix form

xmidmat = [];
ymidmat = [];
cosi = [];
sini = [];

for i = 1:num-1;
    xmidmat = [xmidmat xmid'];
    ymidmat = [ymidmat ymid'];
    cosi = [cosi cos(thetpan') ];
    sini = [sini sin(thetpan')];
end

gamj =  -I1ij/2/pi + I2ij/2/pi;
gamjp1 =  - I2ij/2/pi;


veczero = zeros(1,num-1);

gamjmat = [gamj veczero'];
gamjp1mat = [veczero' gamjp1];

% dphi/dy for imposed vortex
vor_phi_y = (gam_imp/(2*pi))*((xnew(1:end-1)-x_imp)./((xnew(1:end-1)-x_imp).^2+(ynew(1:end-1)-y_imp).^2));

% dphi/dx for imposed vortex
vor_phi_x = (gam_imp/(2*pi))*((ynew(1:end-1)-y_imp)./((xnew(1:end-1)-x_imp).^2+(ynew(1:end-1)-y_imp).^2));

mattot = gamjmat + gamjp1mat;

% set up right hand side

% delta = vinf*h*cos(alp);
% gup = vinf*h*sin(alp);

rhs1 =  vinf*sin(thetpan - alp)   +vor_phi_x.*sin(thetpan)-  vor_phi_y.*cos(thetpan);     % Add vortex equations here ??


        
% add Kutta condition -- last line of matrix and RSH
lastline = zeros(1,num);
lastline(1) = .1;
lastline(end) = .1;
rhs =  [rhs1 0];
 
mat = [mattot ; lastline];

gamans = inv(mat)*rhs';