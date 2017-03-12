% run find_gams_justaf_fix first ... then this computes phi and psi on the
% surface at the midpoints. 

clear xpmat ypmat upan vpan psiloc philoc  psilocpan philocpan philocnew


% xpmat = xmidmat(1:1:end,1); 
% ypmat = ymidmat(1:1:end,1) ;

% 
xpmat = xmidmat(1:1:end,1) -.000001*sini(1:1:end,1);
ypmat = ymidmat(1:1:end,1)  + .000001*cosi(1:1:end,1);

% 
% 
 numx = length(xmid);
 numy = length(ypmat);



for j = 1:1
    for i = 1:numy
        A1 =  -(xpmat(i,j) - xnew(1:end-1)).*cos(thetpan) - (ypmat(i,j) - ynew(1:end-1)).*sin(thetpan);
        B1 =  (xpmat(i,j) - xnew(1:end-1)).^2 + (ypmat(i,j) - ynew(1:end-1)).^2 ;
        E1 =  (xpmat(i,j) - xnew(1:end-1)).*sin(thetpan) - (ypmat(i,j) - ynew(1:end-1)).*cos(thetpan);
       
        temp1 = 2*A1.*SJ + SJ.^2 + B1;
        

           
                ga = atan2(E1.*SJ,(B1+SJ.*A1))./E1;
                 gb = 1/2*log(temp1./B1)  - A1./E1.*atan2(E1.*SJ,(B1+A1.*SJ));
                gc = SJ - A1.*log(temp1./B1)  - atan2(E1.*SJ,(B1+A1.*SJ)) .*(-2*A1.^2 + B1)./E1;
         
                % take care of any NANs
                
                
                
                
%                 
           for iplace = 1:length(E1)
              
                if (abs(E1(iplace)) < 1e-9)
          
                ga(iplace) = 0;
               gb(iplace) = 1/2*log(temp1(iplace)/B1(iplace));
               gc(iplace) = SJ(iplace) - A1(iplace).*log(temp1(iplace)./B1(iplace)) ;
                end
           end
           
           
           % get panel contributions to the velocity field
           
            up = gamans(1:end-1)'./2/pi.* (  (-ypmat(i,j)+ynew(1:end-1)).*(ga) + sin(thetpan).*(gb)  )+ ...
               (gamans(2:end)'-gamans(1:end-1)')./2/pi./SJ.* ( ( -ypmat(i,j) + ynew(1:end-1)).*(gb) + sin(thetpan).*(gc));
           upan(i,j) = -sum(up);
           
           
%             upj(i,:) = 1/2/pi.* (  (-ypmat(i,j)+ynew(1:end-1)).*(ga) + sin(thetpan).*(gb)  )+ ...
%                -(1)./2/pi./SJ.* ( -ypmat(i,j) + ynew(1:end-1)).*(gb) + sin(thetpan).*(gc);
%            upjp1(i,:)= (1)./2/pi./SJ.* ( -ypmat(i,j) + ynew(1:end-1)).*(gb) + sin(thetpan).*(gc);
           
           upj = gamans(1:end-1)'/2/pi.* (  (-ypmat(i,j)+ynew(1:end-1)).*(ga) + sin(thetpan).*(gb)  )+ ...
               (-gamans(1:end-1)')./2/pi./SJ.* ( -ypmat(i,j) + ynew(1:end-1)).*(gb) + sin(thetpan).*(gc);
           upjp1 = (gamans(2:end)')./2/pi./SJ.* ( -ypmat(i,j) + ynew(1:end-1)).*(gb) + sin(thetpan).*(gc);
           
           
           vp = gamans(1:end-1)'./2/pi.* (  (xpmat(i,j) - xnew(1:end-1)).*(ga) - cos(thetpan).*(gb) ) + ...
               (gamans(2:end)'-gamans(1:end-1)') ./2/pi./SJ.*  (  (xpmat(i,j) - xnew(1:end-1)).*(gb) - cos(thetpan).*(gc));
          
           
           vpj =  gamans(1:end-1)'/2/pi.* (  (xpmat(i,j) - xnew(1:end-1)).*(ga) - cos(thetpan).*(gb) ) + ...
               (-gamans(1:end-1)') ./2/pi./SJ.*  (  (xpmat(i,j) - xnew(1:end-1)).*(gb) - cos(thetpan).*(gc));
           vpjp1 = (gamans(2:end)') ./2/pi./SJ.*  (  (xpmat(i,j) - xnew(1:end-1)).*(gb) - cos(thetpan).*(gc));
%           
%              vpj(i,:) =  1/2/pi.* (  (xpmat(i,j) - xnew(1:end-1)).*(ga) - cos(thetpan).*(gb) )  ...
%                -1 ./2/pi./SJ.*  (  (xpmat(i,j) - xnew(1:end-1)).*(gb) - cos(thetpan).*(gc));
%            vpjp1(i,:) = 1 ./2/pi./SJ.*  (  (xpmat(i,j) - xnew(1:end-1)).*(gb) - cos(thetpan).*(gc));
           vpan(i,j) = -sum(vp);
           
           
%          upj*(-sin(thetpan(i))) + vpj*(cos(thetpan(i)))
%          gamj(i,j) 
           
     
          
           % get source, sink and vortices contribution to velocity field
           
           upan(i,j) = upan(i,j) + vinf*cos(alp);
           
           vpan(i,j) = vpan(i,j) + vinf * sin(alp);
           
         
          
    end
end
%     figure(21)
%         hold on;
%     
%         
%         upan.*(-sin(thetpan')) + vpan.*(cos(thetpan'))
% % matA = [upj veczero'];
% % matB = [veczero' upjp1];
% % matC = [vpj veczero'];
% % matD = [veczero' vpjp1];
% %matE = [(matA + matB).*(-sin(thetpan')) + (matC+matD).*cos(thetpan')];
% 
%         %quiver( xpmat, ypmat, upan, vpan,1)  
        
   % Obtain end velocity at trailing edge, in order to determine spacing for
   % first shed vortex
 
   Vx = upan(1,1)+ upan(end);
   %Vy = vpan(1,1)+vpan(end);
   V_end = Vx/2;
   
% Say spacing in time dt = 0.1 for now, vortex will be sheded at 
dt = 0.01;
x_fixed(1) = xnew(end) +0.504*dt;
%x_fixed(1) = 0.55;
%x_fixed(1) = 0.5+ xl(end) - xl(end-1);    
    
        
    cp = 1 -((upan.^2 + vpan.^2)); 
    figure(5)
    plot(xpmat,-cp);
    hold on
    cp2 = 1 - ((upan.^2 + vpan.^2));
    plot(xpmat, -cp2,'r');
    
    
%    % compute the lift coefficient -- in calculations chord is 1 (-.5 to .5)
   lift_coef = sum(cp'.*SJ.*cos(thetpan))/vinf^2;
% 
%     % when computing velocity field on body - use this to check normal
%     % condition
%     
%         for i = 1:length(xpmat)
%             checker(i) = upan(i,1)*(-sin(thetpan(i))) + vpan(i,1)*cos(thetpan(i)) ;
%         end;
%         
%         [maxxnew locxmax] = max(xnew);
% [minxnew locxmin] = min(xnew);
%         %Build potential field due to panels
%         
%         % Loop over points where you want the value of the potential
%         for j = 1:1
%             for i = 1:numy
%              
%                 
%                 a = ypmat(i,j) - ynew(1:end-1);
%                 b = -sin(thetpan);
%                 c = xpmat(i,j) - xnew(1:end-1);
%                 d = -cos(thetpan);
%                 t = (a.*d-b.*c);
%                 m = b./d + t./(d.* ( d.*SJ + c) );
%                 n = (b.^2.*SJ + b.*a + d.^2.*SJ + d.*c)./t;
%                 p = t./(d.*SJ + c);
%                 q = 1/2*a.^2.*b.^2 - 1/2*d.^2.*a.^2 - 3/2.*b.^2.*c.^2 + 2*a.*b.*c.*d - 1/2.*b.^4.*c.^2./d.^2;
%                 s = d.^2.*c.*a - d.*c.^2.*b + d.*a.^2.*b - a.*b.^2.*c;
%                 fac1 = (d.*SJ + c) ./d;
%                 
%                 f1 = b.^2.*SJ.^2 + d.^2.*SJ.^2 + 2*a.*b.*SJ + 2*c.*d.*SJ + a.^2 + c.^2;
%                 g1 = (SJ + a.*b + c.*d);
%                 
%                 
%                A1 =  -(xpmat(i,j) - xnew(1:end-1)).*cos(thetpan) - (ypmat(i,j) - ynew(1:end-1)).*sin(thetpan);
%         B1 =  (xpmat(i,j) - xnew(1:end-1)).^2 + (ypmat(i,j) - ynew(1:end-1)).^2 ;
%         E1 =  (xpmat(i,j) - xnew(1:end-1)).*sin(thetpan) - (ypmat(i,j) - ynew(1:end-1)).*cos(thetpan);
%        
%         temp1 = 2*A1.*SJ + SJ.^2 + B1;
%         
%                  ga = (atan2(E1.*SJ,(B1+SJ.*A1))); 
%                 gb = log(temp1./B1);
%                 gc =  (atan2(a + b.*SJ, c + d.*SJ));
%               
% %                 while ( max((gc(2:end)-gc(1:end-1))) > 3  | min((gc(2:end)-gc(1:end-1))) < -3) 
% %                     [val loc] =      max((gc(2:end)-gc(1:end-1))) ;
% %                     while val > 3
% %                         if  abs(gc(loc)) > abs(gc(loc+1))
% %                            gc(loc) = gc(loc) + pi;
% %                        else
% %                           
% %                             gc(loc+1) = gc(loc+1) - pi;
% %                         end
% %                             [val loc] =      max((gc(2:end)-gc(1:end-1))) ;
% %                     end
% %                     [val loc] = min((gc(2:end)-gc(1:end-1))) ;
% %                     while val < -3
% %                         if abs(gc(loc)) > abs(gc(loc+1))
% %                             gc(loc)  = gc(loc) -pi;
% %                         else
% %                         gc(loc+1) = gc(loc+1) + pi;
% %                         end
% %                         [val loc] = min((gc(2:end)-gc(1:end-1))) ;
% %                     end
% %                 end
% %                             
% %                 
% %                 while ( max((ga(2:end)-ga(1:end-1))) > 3  | min((ga(2:end)-ga(1:end-1))) < -3) 
% %                     [val loc] =      max((ga(2:end)-ga(1:end-1))) ;
% %                     while val > 3
% %                        if  abs(ga(loc)) > abs(ga(loc+1))
% %                            ga(loc) = ga(loc) + pi;
% %                        else
% %                           
% %                             ga(loc+1) = ga(loc+1) - pi;
% %                        end
% %                             [val loc] =      max((ga(2:end)-ga(1:end-1))) ;
% %                     end
% %                     [val loc] = min((ga(2:end)-ga(1:end-1))) ;
% %                     while val < -3
% %                         if abs(ga(loc)) > abs(ga(loc+1))
% %                             ga(loc)  = ga(loc) -pi;
% %                         else
% %                             
% %                         ga(loc+1) = ga(loc+1) + pi;
% %                         end
% %                         [val loc] = min((ga(2:end)-ga(1:end-1))) ;
% %                      
% %                     end
% %                 end
%                 
%                  
% %              
%                 
%                 
%                 %  y - ymin = (ymax - ymin)/(xmax-xmin)*(x - xmin)
%                         yonchord = (ynew(locxmax) - ynew(locxmin))/(xnew(locxmax)-xnew(locxmin))*(xpmat(i,j) - xnew(locxmin)) + ynew(locxmin);
%                        
%                 if ( xpmat(i,j) > minxnew & xpmat(i,j) < maxxnew & ypmat(i,j) < yonchord |  i < locxmin) 
%                     
%                     
%                     for kk = 1:length(ynew)-1
%                     if (a(kk) + b(kk).*SJ(kk) < 0) 
%                             gc(kk) = gc(kk) + 2*pi;
%                     end 
%                     if (a(kk) + b(kk).*SJ(kk) > 0 & c(kk) + d(kk).*SJ(kk) > 0) 
%                         gc(kk) = gc(kk) + 2*pi;
%                     end
%                     end
%                     
%                 else   
%                     
%                 for kk = 1:length(ynew)-1
%                     if (a(kk) + b(kk).*SJ(kk) < 0) 
%                         if ( c(kk) + d(kk).*SJ(kk)  < 0 )
%                             gc(kk) = gc(kk) + 2*pi;
%                         end
%                        
%                     end
%                     
%                 end
%                 end
%                   
%                 
%                 ztop = (1/2.*d.^2.*a.^2 - a.*b.*c.*d - c.*a.*d.*p + c.^2.*b.*p + 1/2*b.^2 .* c.^2);
%                 wtop = -a.*b.*c.*d + 1/2.*b.^2.*c.^2 + 1/2.*d.*a.^2.*d;
%                 
%                 fac2 = (-a.^3.*b-a.^2.*c.*d-a.*b.*c.^2-c.^3.*d+a.^3.*b.^3-a.^3.*b.*d.^2+5*a.^2.*b.^2.*c.*d - a.^2.*c.*d.^3 - a.*b.^3.*c.^2 ...
%                     + 5.*a.*b.*c.^2.*d.^2 - b.^2.*c.^3.*d + c.^3.*d.^3);
%       
%                 fac3 = (a.^2 + c.^2 - a.^2.*b.^2 - 2*a.*b.*c.*d - c.^2.*d.^2);
%                
%                 for jj = 1:length(ztop)
%                 if (abs(ztop(jj)) < 1e-12)
%                         z(jj) = 0;
%                 else
%                         z(jj) = ztop(jj)/d(jj)^2/p(jj)^2;
%                 end
%                 
%                 if ( abs(wtop(jj)) < 1e-12  ) 
%                        w(jj) = 0;
%                 else
%                      w(jj) = wtop(jj)/d(jj)/p(jj);
%                 end
%                     
%              %   if     (  abs(fac2(jj)) < 1e-8)
%                 if (abs(t(jj)) < 1e-12)
%                     
%                     f3(jj) = 0;
%                     f4(jj) = 0;
%                    
%                 else 
%                    
%                     f3(jj) = fac2(jj).*atan(g1(jj)/t(jj))./t(jj);
%                      f4(jj) = 2*atan(g1(jj)/t(jj)).*fac3(jj)./t(jj);
%                 end
%                end
%                 
%           
%                 inter1_SJ = atan(m).*fac1  + t.*log(m.^2+1)./2  ...
%                      - atan(n).*b.*t./d;
%                 inter2_SJ =  + atan(n).*q  - log(m.^2+1).*s./2 +atan(m).*z + w;
%                 
%                 psiint1_SJ = log(f1).*(SJ + a.*b + c.*d) + f4 -2*SJ  ;
%                 psiint2_SJ = 1/2*log(f1).*(SJ.^2 - a.^2.*b.^2+a.^2.*d.^2-4*a.*b.*c.*d+b.^2.*c.^2-c.^2.*d.^2) + f3 ...
%                      - 1/2*SJ.^2 + a.*b.*SJ + c.*d.*SJ;
%                 
%           
%                 m = b./d + t./(d.* c );
%                 n = ( b.*a + d.*c)./t;
%                 p = t./c;
%                 fac1 = ( c) ./d;
%                  
%                 f1 =  a.^2 + c.^2;
%                 g1 = (a.*b + c.*d);
%                 
%                 
%                 ztop = (1/2.*d.^2.*a.^2 - a.*b.*c.*d - c.*a.*d.*p + c.^2.*b.*p + 1/2*b.^2 .* c.^2);
%                  wtop = -a.*b.*c.*d + 1/2.*b.^2.*c.^2 + 1/2.*d.*a.^2.*d;
%       
%                  
%                for jj = 1:length(ztop)
%                   
%                    
%                 if (abs(ztop(jj)) < 1e-12)
%                         z(jj) = 0;
%                 else
%                         z(jj) = ztop(jj)/d(jj)^2/p(jj)^2;
%                 end
%                  if (abs(wtop(jj)) < 1e-12) 
%                        w(jj) = 0;
%                 else
%                      w(jj) = wtop(jj)/d(jj)/p(jj);
%                  end
%                   
%                 if (abs(t(jj)) < 1e-12)
%                     
%                     f3(jj) = 0;
%                     f4(jj) = 0;
%                    
%                 else 
%                    
%                     f3(jj) = fac2(jj)*atan(g1(jj)/t(jj))./t(jj);
%                      f4(jj) = 2*atan(g1(jj)/t(jj)).*fac3(jj)./t(jj);
%                 
%                end
%                end
%                
%                %%%%
%       % new phi way - totally integrated
%       %%%%
%                  phiI = gc.*SJ  + E1/2.*gb   -A1.* ga;
%                 i;
%                 j;
%                 figure(9)
%                 hold off
%                 plot(phiI)
%                 hold on
%                 axis([ 0 161  -.2 .2])
%                 
%                
%                phiII_overS = E1/2 -ga.*(-2*A1.^2 + B1)./SJ/2 - E1.*A1.*gb./SJ/2 + SJ/2.*gc;
%                plot(phiII_overS)
%               % pause
%                 
%                 ansphipan = gamans(1:end-1)'.*(phiI)/2/pi  ...
%                    +    ( gamans(2:end)' - gamans(1:end-1)'  ).*(phiII_overS)/2/pi  ;
%                 philocpannew(i,j) = - sum(real(ansphipan));
%                 philocnew(i,j) = -sum(real(ansphipan))  + vinf*cos(alp)*xpmat(i,j) + vinf*sin(alp)*ypmat(i,j);
%                 %%%%%%%%%%%%%
%                 
%                 inter1_0 =  atan(m).*fac1  + t.*log(m.^2+1)./2  ...
%                     - atan(n).*b.*t./d;
%                 inter2_0 =  + atan(n).*q  - log(m.^2+1).*s./2 +atan(m).*z + w;
%                 
%                  psiint1_0 = log(f1).*(a.*b + c.*d) + f4  ;
%                 psiint2_0 = 1/2*log(f1).*(- a.^2.*b.^2+a.^2.*d.^2-4*a.*b.*c.*d+b.^2.*c.^2-c.^2.*d.^2) + f3 ;
%                      
%                 
%                 
%                 ansint1 = inter1_SJ - inter1_0 - t.*log(c./(d.*SJ+c));
%                 ansint2 = inter2_SJ - inter2_0 + s.*log(c./(d.*SJ + c));
%                 
%                 psiansint1 = psiint1_SJ - psiint1_0;
%                 psiansint2 = psiint2_SJ - psiint2_0;
%                 
%                 anspsipan = gamans(1:end-1)'.*(psiansint1)/4/pi  ...
%                    +    ( gamans(2:end)' - gamans(1:end-1)'  )./SJ.*(psiansint2)/4/pi  ;
%                 
%                 
%                 
%                 anspan = gamans(1:end-1)'.*(ansint1)/2/pi  ...
%                    +    ( gamans(2:end)' - gamans(1:end-1)'  )./SJ.*(ansint2)/2/pi  ;
%                 
%                
%                philocpan(i,j) = -sum(real(anspan));
%                 philoc(i,j) = -sum(real(anspan)) + vinf*cos(alp)*xpmat(i,j) + vinf*sin(alp)*ypmat(i,j);
%                 
%                 
%                 psilocpan(i,j) = sum(anspsipan);
%                 psiloc(i,j) = sum(anspsipan) + vinf*cos(alp)*ypmat(i,j) - vinf*sin(alp)*xpmat(i,j);
%                 
% 
%             end
%         end
%         
%         figure(20)
%         plot(xpmat,(philocnew))
%         
%          figure(22)
%         plot(xpmat,psiloc)
%         hold on
%         plot(xpmat,philocnew,'r');
        
        
       
        
        
%         figure(23)
%         pcolor(xpmat,ypmat,psiloc);
%         hold on
%         plot(xnew,ynew);
%         
%         figure(24);
%         contour(xpmat,ypmat,psiloc);
%         hold on
%         contour(xpmat,ypmat,philoc);
%         plot(xnew,ynew)
%         