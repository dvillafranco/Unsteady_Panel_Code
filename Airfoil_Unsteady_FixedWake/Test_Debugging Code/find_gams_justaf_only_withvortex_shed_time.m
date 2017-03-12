% find_gams_justaf_only_withvortex


for n = 1:3
% set up matrix to find unknown gamma values 
% use casc_tranform.m to get the endpoints, midpoints, and panel angles and
%   h which is the cascade gap

% Define Parameters for the Vortex
gam_imp =0.5;
x_imp = -2.5;
y_imp = 0.03;

% Define Parameters for Shed Vorted
%x_fixed = 0.6;
y_fixed = 0.03;
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

x_vor = x_fixed(n);

% dphi/dy for imposed vortex
vor_phi_y = (-gam_imp/(2*pi))*((xnew(1:end-1)-x_imp)./((xnew(1:end-1)-x_imp).^2+(ynew(1:end-1)-y_imp).^2));

% dphi/dx for imposed vortex
vor_phi_x = (gam_imp/(2*pi))*((ynew(1:end-1)-y_imp)./((xnew(1:end-1)-x_imp).^2+(ynew(1:end-1)-y_imp).^2));

% dphi/dy for shedded vortex
vor_phi_yshed = -1*((-1/(2*pi))*((xnew(1:end-1)-x_vor)./((xnew(1:end-1)-x_vor).^2+(ynew(1:end-1)-y_fixed).^2))).*cos(thetpan);

% dphi/dx for shedded vortex
vor_phi_xshed = ((1/(2*pi))*((ynew(1:end-1)-y_fixed)./((xnew(1:end-1)-x_vor).^2+(ynew(1:end-1)-y_fixed).^2))).*sin(thetpan);



vormat(:,n) = (vor_phi_yshed-vor_phi_xshed)';


mattot = [(gamjmat + gamjp1mat), vormat];

% set up right hand side

% delta = vinf*h*cos(alp);
% gup = vinf*h*sin(alp);

if n == 1
    rhs1 =  vinf*sin(thetpan - alp)   +vor_phi_x.*sin(thetpan)-  vor_phi_y.*cos(thetpan);%+vor_phi_xshed  - vor_phi_yshed;     % Add vortex equations here ??
else
    rhs1 = vinf*sin(thetpan - alp)   +vor_phi_x.*sin(thetpan)-  vor_phi_y.*cos(thetpan) - vormat(:,n-1)';
end

steady_mat = gamans_steady.*ones(1,n);
        
% add Kutta condition -- last line of matrix and RSH
lastline = zeros(1,length(mattot));
lastline(1) = .1;
lastline(end) = .1;
rhs =  [rhs1 steady_mat 0  ];


% Kelvin's Theory to keep circulation sum zero
kelvin = ones(1,length(mattot));

LL = repmat(kelvin,n,1);

 
mat = [mattot ;LL; lastline];

gamans = inv(mat)*rhs';

%gamans_matrix(161,n) = gamans;

x_fixed(n+1) = x_fixed(n)*2;

end


