% set up matrix to find unknown gamma values 
% use casc_tranform.m to get the endpoints, midpoints, and panel angles and
%   h which is the cascade gap


num = length(xnew);


delta = vinf*h*cos(alp);
gup = vinf*h*sin(alp);

% delta = vinf*cos(alp);
%  gup = vinf*sin(alp);

% set up constants that appear in the integrals

SJmat = [];
for i = 1:num-1;
        SJmat = [SJmat;  SJ];
end;

for i = 1:num-1
        A(i,:) = -(xmid(i) - xnew(1:end-1)).*cos(thetpan) - (ymid(i) - ynew(1:end-1)).*sin(thetpan);
        B(i,:) = (xmid(i) - xnew(1:end-1)).^2 + (ymid(i) - ynew(1:end-1)).^2;
        D(i,:) = cos(thetpan(i) - thetpan);
        C(i,:) = -sin(thetpan(i) - thetpan);
        E(i,:) = (xmid(i) - xnew(1:end-1)).*sin(thetpan) - (ymid(i) - ynew(1:end-1)).*cos(thetpan);
        R(i,:) = (xmid(i) - xnew(1:end-1))*cos(thetpan(i)) + (ymid(i) - ynew(1:end-1)).*sin(thetpan(i));
        T11(i,:) = -(xmid(i) - xnew(1:end-1))/2.*cos(thetpan(i) - 2*thetpan) ...
            + (ymid(i) - ynew(1:end-1))/2.*sin(thetpan(i)-2*thetpan) ;
        T2(i,:) = - (  (xmid(i) - xnew(1:end-1)).*sin(thetpan)-(ymid(i) - ynew(1:end-1)).*cos(thetpan) ).*cos(thetpan(i)-thetpan) ...
            +   (   (xmid(i) - xnew(1:end-1)).*cos(thetpan)-(ymid(i) - ynew(1:end-1)).*sin(thetpan)).*sin(thetpan(i)-thetpan);
        
        Gd1(i,:) = SJmat(i,:)/2*( ymid(i)/ ( (xmid(i) -1)^2 + ymid(i)^2) )*(-sin(thetpan(i))) ;
        Gd2(i,:) = SJmat(i,:)/2*( (xmid(i)-1)/( (xmid(i) -1)^2 + ymid(i)^2)) *cos(thetpan(i));
        
end;

% create SJ matrix (repeat in each column)



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

gamj =  -I1ij/2/pi + I2ij/2/pi  - Gd1/2/pi + Gd2/2/pi;
gamjp1 =  - I2ij/2/pi - Gd1/2/pi  + Gd2/2/pi;


veczero = zeros(1,num-1);

gamjmat = [gamj veczero'];
gamjp1mat = [veczero' gamjp1];



mattot = gamjmat + gamjp1mat;

% set up right hand side

bot3 = (xmid - 1).^2 + ymid.^2;
bot4 = (xmid + 1).^2 + ymid.^2;

rhs1 = delta/2/pi*(  ((xmid+1)./bot4 - (xmid-1)./bot3).*sin(thetpan) + (ymid./bot4 - ymid./bot3).*(-cos(thetpan) )  ) ...
            -  gup /2/pi.*(  (ymid./bot4 - ymid./bot3).*sin(thetpan) + ( - (xmid+1)./bot4  + (xmid-1)./bot3) .*(-cos(thetpan)) );        

        
% add Kutta condition -- last line of matrix and RSH
lastline = zeros(1,num);
lastline(1) = .1;
lastline(end) = .1;
rhs =  [rhs1 0];

mat = [mattot ; lastline];

gamans = inv(mat)*rhs'


for j = 1:length(SJ)
    integral(j) = (gamans(j) + gamans(j+1))*SJ(j)/2;
end

gamans_steady = sum(integral);


        
            
        