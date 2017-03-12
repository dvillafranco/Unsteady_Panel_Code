% find_gams_justaf_only_withvortex
ntz=1;
for nx = 1:1200
% set up matrix to find unknown gamma values 
% use casc_tranform.m to get the endpoints, midpoints, and panel angles and
%   h which is the cascade gap

% Define Parameters for the Vortex
gam_imp =0.01;
x_imp(1) = -9.5;
y_imp(1) = 0.1;

% Define Parameters for Shed Vorted
%x_fixed = 0.6;
y_fixed = 0.0;
%gam_shed = 0.2;

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
vor_phi_y = (gam_imp/(1))*((xnew(1:end-1)-x_imp(nx))./((xnew(1:end-1)-x_imp(nx)).^2+(ynew(1:end-1)-y_imp(1)).^2));

% dphi/dx for imposed vortex
vor_phi_x = (-gam_imp/(1))*((ynew(1:end-1)-y_imp(1))./((xnew(1:end-1)-x_imp(nx)).^2+(ynew(1:end-1)-y_imp(1)).^2));


% dphi/dy for shedded vortex
vor_phi_yshed = ((1/(1))*((xnew(1:end-1)-x_fixed(1))./((xnew(1:end-1)-x_fixed(1)).^2+(ynew(1:end-1)-y_fixed).^2))).*cos(thetpan);

% dphi/dx for shedded vortex
vor_phi_xshed = ((1/(1))*((ynew(1:end-1)-y_fixed)./((xnew(1:end-1)-x_fixed(1)).^2+(ynew(1:end-1)-y_fixed).^2))).*sin(thetpan);


%mattot = [(gamjmat + gamjp1mat), (vor_phi_yshed+vor_phi_xshed)'];

mattot = [(gamjmat + gamjp1mat), (-vor_phi_yshed - vor_phi_xshed)'];

%mattot = [(gamjmat + gamjp1mat), zeros(size(thetpan))'];


% set up right hand side

% delta = vinf*h*cos(alp);
% gup = vinf*h*sin(alp);



if nx == 1
    rhs1 =  vinf*sin(thetpan - alp)   -vor_phi_x.*sin(thetpan)+  vor_phi_y.*cos(thetpan);%+vor_phi_xshed  - vor_phi_yshed;     % Add vortex equations here ??
else
    rhs1 = vinf*sin(thetpan - alp)   -vor_phi_x.*sin(thetpan)+  vor_phi_y.*cos(thetpan) + sum(gam_wake,2)'; %vormat(:,n-1)';
end


%rhs1 =  vinf*sin(thetpan - alp)   +vor_phi_x.*sin(thetpan)-  vor_phi_y.*cos(thetpan) +vor_phi_xshed  - vor_phi_yshed;     % Add vortex equations here ??


        
% add Kutta condition -- last line of matrix and RSH
lastline = zeros(1,num+1);
lastline(1) = 1;
lastline(end-1) = 1;

if nx == 1
    rhs =  [rhs1  0 gamans_steady];%/(2*pi)];
else 
    rhs =  [rhs1 0 (gamans_steady-sum(gam_dim))];
end

% Kelvin's Theory to keep circulation sum zero
kelvin = ones(1,num+1);
 
mat = [mattot ; lastline; kelvin];

gamans = inv(mat)*rhs';

gamans_MAT(:,nx) = gamans;

if nx == 1 
    gam_dim = gamans(end);
else 
   gam_dim = [gamans(end), gam_dim];
end 

%x_fixed(n+1) = x_fixed(1)*(n+1);

% for aa = length(gam_dim)
%     gam_wake(:,aa) = (1*gam_dim(aa).*((1/(1))*((xnew(1:end-1)-x_fixed(aa))./((xnew(1:end-1)-x_fixed(aa)).^2+(ynew(1:end-1)-y_fixed).^2))).*cos(thetpan)+...
%     (1*gam_dim(aa)*((1/(1))*((ynew(1:end-1)-y_fixed)./((xnew(1:end-1)-x_fixed(aa)).^2+(ynew(1:end-1)-y_fixed).^2))).*sin(thetpan)));
% end

psionsurface;

%psiinfield_vortex;


ntz = ntz+1

end
