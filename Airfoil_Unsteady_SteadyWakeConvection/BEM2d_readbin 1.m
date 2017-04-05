function D = BEM2d_readbin(filetag)
%D = BEM2d_readbin('filetag')  reads binary data from BEM.filetag.bin
%
% Matlab post processing code for BEM2d2.f90
%	AUTHOR:  Trevor Wood, Boston University
%	DATE: 	c. 2000
%
%The output structure D contains the following data:
%
%   t       time vector
%   c       airfoil chord length
%   l       airfoil panel lengths
%   n       outward normal vectors of airfoil panels
%   Xb      airfoil nodes           [xn,yn]
%   Xc      collocation points      [xc,yc]
%   Xw      wake nodes at t(end)    [xw,yw]
%   Pb      perturbation velocity potential, phi, at (Xc_i,t_j)
%   Pbn     phi at Xb using spline inter/extrapolation
%   Pw      potential jump of wake panels
%   phidot  time derivative of phi (using 2nd order forw/cent/back)
%   vt      tangential flow velocity at Xc
%   Cp      Coefficient of pressure
%   L       lift coefficient scaled by (rho U^2 c/2)
%   M       lift moment coef about +x_2 axis scaled by (rho U^2 c^2/2)
%
%

if nargin<1 
   filetag = ''; 
else
   filetag = [filetag,'.'];
end

file = ['BEM_data.',filetag,'bin'];
fid = fopen(file,'r');
   D(1).t = readbin(fid,1);
   Ne     = readbin(fid);

   for i = 1:Ne
      D(i).Xb = readbin(fid,2);
      D(i).Xc = (D(i).Xb(2:end,:) + D(i).Xb(1:end-1,:))/2;
   end
   Xc = readbin(fid,2);

   for j = 1:length(D(1).t)
   for i = 1:Ne
      D(i).Pb(:,j) = readbin(fid,1);
      D(i).Xw{j}   = readbin(fid,2);
      D(i).Pw{j}   = readbin(fid,1);
   end
   end
fclose(fid);

D = getVT(D);

%D = getdphi1(D);
D = getdphi2(D);

for j=1:length(D)
   D(j).Cp = 1 - 2*D(j).phidot - D(j).vt.^2;
%   D(j).Cp = 1 - D(j).vt.^2;
%   D(j).Cp =  - 2*D(j).phidot ;
for i=1:length(D(1).t)
   D(j).L(i) = sum(-D(j).Cp(:,i).*D(j).n(:,2).*D(j).l) / D(j).c;
   D(j).M(i) = sum( D(j).Cp(:,i).*D(j).n(:,2).*D(j).l.*D(j).Xc(:,1)) / D(j).c^2;
end
end


function D = getVT(D)
   for i = 1:length(D)
      D(i).c = max(sqrt((D(i).Xb(2:end,1)-D(i).Xb(1,1)).^2+(D(i).Xb(2:end,2)-D(i).Xb(1,2)).^2));
      D(i).l = sqrt(sum((D(i).Xb(2:end,:)-D(i).Xb(1:end-1,:)).^2,2));
      t      = (D(i).Xb(2:end,:)-D(i).Xb(1:end-1,:))./[D(i).l D(i).l];
      D(i).n = [-t(:,2) t(:,1)];

      sn = [0 D(i).l'];
      for j=1:length(D(i).l)
         s(j)  = sum(sn(1:j));
         sc(j) = s(j) + 0.5*sn(j+1);
      end
      s(length(D(i).l)+1) = sum(sn);

      D(i).Pbn = spline(sc,D(i).Pb',s)';
      D(i).vt  = (D(i).Pbn(2:end,:)-D(i).Pbn(1:end-1,:))./(D(i).l*ones(1,length(D(1).t)));
      D(i).vt  = D(i).vt + (t(:,1)*ones(1,length(D(1).t)));
   end

function D = getdphi1(D)
   dt = diff(D(1).t);
   if(length(dt) < 2)
      for i = 1:length(D)
         D(i).phidot = 0*D(i).Pb;
      end
   end

   if any(abs(dt-dt(1)) > 1e-3*dt(1)) error 'Data not equally spaced'; end;
   dt = dt(1);

   for i = 1:length(D)
      D(i).phidot = 0*D(i).Pb;
      D(i).phidot(:,2:end) = ( D(i).Pb(:,2:end) - D(i).Pb(:,1:end-1) ) / dt;
      D(i).phidot(:,  1  ) = D(i).Pb(:,1) / dt;
   end

function D = getdphi2(D)
   dt = diff(D(1).t);
   if(length(dt) < 2)
      for i = 1:length(D)
         D(i).phidot = 0*D(i).Pb;
      end
   end

   if any(abs(dt-dt(1)) > 1e-3*dt(1)) error 'Data not equally spaced'; end;
   dt = dt(1);
   
   for i = 1:length(D)
      D(i).phidot = 0*D(i).Pb;
   D(i).phidot(:,2:end-1) = ( D(i).Pb(:,3:end) -   D(i).Pb(:,1:end-2) );
      D(i).phidot(:,   1   ) = (-D(i).Pb(:,  3  ) + 4*D(i).Pb(:,   2   ) - 3*D(i).Pb(:, 1 ) );
      D(i).phidot(:,  end  ) = ( D(i).Pb(:,end-2) - 4*D(i).Pb(:, end-1 ) + 3*D(i).Pb(:,end) );
      D(i).phidot            = D(i).phidot / (2*dt);
   end
