% set up matrix to find unknown gamma values 
% use casc_tranform.m to get the endpoints, midpoints, and panel angles and
%   h which is the cascade gap


num = length(xnew);

% vinf = 1;  % freestream velocity
% alp = 0*pi/180;   % freestream angle of attack

% set up constants that appear in the integrals

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

mattot = gamjmat + gamjp1mat;

% set up right hand side

% delta = vinf*h*cos(alp);
% gup = vinf*h*sin(alp);

rhs1 =  vinf*sin(thetpan - alp);


        
% add Kutta condition -- last line of matrix and RSH
lastline = zeros(1,num);
lastline(1) = .1;
lastline(end) = .1;
rhs =  [rhs1 0];
 
mat = [mattot ; lastline];

gamans = inv(mat)*rhs';

%gamans_steady = sum(gamans);

for j = 1:length(SJ)
    integral(j) = (gamans(j) + gamans(j+1))*SJ(j)/2;
end

gamans_steady = sum(integral);





