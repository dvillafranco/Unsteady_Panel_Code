% run find_gams_justaf_fix first ... then this computes phi and psi on the
% surface at the midpoints. 

clear xpmat ypmat upan vpan psiloc philocnew  psilocpan  psilocpan philocpan philocnew

xpmat = xmidmat(1:1:end,1)-.000001*sini(1:1:end,1);
 ypmat = ymidmat(1:1:end,1)+ .000001*cosi(1:1:end,1);

% 
% 
 numx = length(xmid);
 numy = length(ymid);
bot33 = (xpmat - 1).^2 + ypmat.^2;
           bot44 = (xpmat + 1).^2 + ypmat.^2;


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
           for iplace = 1:length(E1)
              
                if (abs(E1(iplace)) < 1e-9)
          
                ga(iplace) = 0;
               gb(iplace) = 1/2*log(temp1(iplace)/B1(iplace));
               gc(iplace) = SJ(iplace) - A1(iplace).*log(temp1(iplace)./B1(iplace)) ;
                end
           end
           
           % get panel contributions to the velocity field
           
            up = gamans(1:end-2)'./2/pi.* (  (-ypmat(i,j)+ynew(1:end-1)).*(ga) + sin(thetpan).*(gb)  )+ ...
               (gamans(2:end-1)'-gamans(1:end-2)')./2/pi./SJ.* ( ( -ypmat(i,j) + ynew(1:end-1)).*(gb) + sin(thetpan).*(gc));
           upan(i,j) = -sum(up);
           
           upj = gamans(1:end-2)'/2/pi.* (  (-ypmat(i,j)+ynew(1:end-1)).*(ga) + sin(thetpan).*(gb)  )+ ...
               (-gamans(1:end-2)')./2/pi./SJ.* ( -ypmat(i,j) + ynew(1:end-1)).*(gb) + sin(thetpan).*(gc);
           upjp1 = (gamans(2:end-1)')./2/pi./SJ.* ( -ypmat(i,j) + ynew(1:end-1)).*(gb) + sin(thetpan).*(gc);
           
           
           vp = gamans(1:end-2)'./2/pi.* (  (xpmat(i,j) - xnew(1:end-1)).*(ga) - cos(thetpan).*(gb) ) + ...
               (gamans(2:end-1)'-gamans(1:end-2)') ./2/pi./SJ.*  (  (xpmat(i,j) - xnew(1:end-1)).*(gb) - cos(thetpan).*(gc));
          
           
           vpj =  gamans(1:end-2)'/2/pi.* (  (xpmat(i,j) - xnew(1:end-1)).*(ga) - cos(thetpan).*(gb) ) + ...
               (-gamans(1:end-2)') ./2/pi./SJ.*  (  (xpmat(i,j) - xnew(1:end-1)).*(gb) - cos(thetpan).*(gc));
           vpjp1 = (gamans(2:end-1)') ./2/pi./SJ.*  (  (xpmat(i,j) - xnew(1:end-1)).*(gb) - cos(thetpan).*(gc));
          
           
           vpan(i,j) = -sum(vp);
           
%          upj*(-sin(thetpan(i))) + vpj*(cos(thetpan(i)));
%          gamj(i,j) ;
           
     
          gamdown = sum( gamans(2:end-1)'.*SJ/2 + gamans(1:end-2)'.*SJ/2) - vinf*h*sin(alp);
           % get source, sink and vortices contribution to velocity field
           
           upan(i,j) = upan(i,j)  + delta/2/pi*  ((xpmat(i,j)+1)./bot44(i,j) - (xpmat(i,j)-1)./bot33(i,j) ) -  gup /2/pi*(ypmat(i,j)./bot44(i,j)) - gamdown/2/pi*(ypmat(i,j)./bot33(i,j))+...
               (gam_imp/(1))*((ypmat(i,j)-y_imp(nx))./((xpmat(i,j)-x_imp(1,nx)).^2+(ypmat(i,j)-y_imp(nx)).^2))+...
               sum((gam_dim./(1)).*((ypmat(i,j)-y_fixed)./((xpmat(i,j)-x_fixed).^2+(ypmat(i,j)-y_fixed).^2)));
           
           
           vpan(i,j) = vpan(i,j) + delta/2/pi*(  ypmat(i,j)./bot44(i,j) - ypmat(i,j)./bot33(i,j)) + gup/2/pi*( (xpmat(i,j)+1)./bot44(i,j) ) + gamdown/2/pi*((xpmat(i,j)-1)./bot33(i,j))+...
               -(gam_imp/(1))*((xpmat(i,j)-x_imp(nx))./((xpmat(i,j)-x_imp(nx)).^2+(ypmat(i,j)-y_imp(nx)).^2))...
               +sum((-gam_dim./(1)).*((xpmat(i,j)-x_fixed)./((xpmat(i,j)-x_fixed).^2+(ypmat(i,j)-y_fixed).^2)));
          
    end
end
        figure(21)
        hold on
        quiver( xpmat, ypmat, upan, vpan,1)
        axis equal
        
%         % Get position for first wake to be shed
%         Vx = upan(1,1) + upan(end);
%         V_end = Vx/2;
%         dt = 0.01;
%         x_fixed(1) = 1 + V_end*dt;
%         %
        
        convect_cascade;
        convect_imposedvor_cascade;
        
for aa = 1:length(gam_dim)
    gam_wake(:,aa) = (1*gam_dim(aa).*((1/(1))*((xnew(1:end-1)-x_fixed(aa+1))./((xnew(1:end-1)-x_fixed(aa+1)).^2+(ynew(1:end-1)-y_fixed(aa+1)).^2))).*cos(thetpan)+...
    (1*gam_dim(aa)*((1/(1))*((ynew(1:end-1)-y_fixed(aa+1))./((xnew(1:end-1)-x_fixed(aa+1)).^2+(ynew(1:end-1)-y_fixed(aa+1)).^2))).*sin(thetpan)));
end

        z2 = xpmat + 1i*ypmat;
        z1 = atanh(z2)*h/pi;
        %z1 = xpinaf + 1i*ypinaf;
        dz2dz1 = pi/h*(1 - (tanh(pi*z1/h)).^2 );
        
        compv2 = upan - 1i*vpan;
        compv1 = compv2.*dz2dz1;
        
        u1 = real(compv1);
        v1 = -imag(compv1);
        
        
        
        
    cpcasc = 1 -(u1.^2 + v1.^2); 
    figure(5)
    plot(real(z1),-cpcasc);
    %axis([-.5 .5 -2 4]);
    title('-Cp')
    
    figure(6)
    plot( real(z1)*cos(chi) + imag(z1)*sin(chi)   ,1-cpcasc)
    % axis([ -.5  .5  0 2.8])
    %set(gca,'XTick',[-.5:0.2:.5])
    %set(gca,'YTick',[0:0.4:2.8])
     grid on
     title('1-Cp')
     
      lift_coef(ntz) = sum(cpcasc'.*SJ.*cos(thetpan))/vinf^2;
     
     
%     % when computing velocity field on body - use this to check normal
%     % condition
%     
%         for i = 1:length(xmid)
%             checker(i) = upan(i,1)*(-sin(thetpan(i))) + vpan(i,1)*cos(thetpan(i)) ;
%         end;
%         
%         
% 
% lastpoints =  round(length(ynew)/2); 
%  [maxxnew locxmax] = max(xnew);
% [minxnew locxmin] = min(xnew);
%                 %Build potential field due to panels
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
%                 
%                 m = b./d + t./(d.* ( d.*SJ + c) );
%                 n = (b.^2.*SJ + b.*a + d.^2.*SJ + d.*c)./t;
%                 p = t./(d.*SJ + c);
%                 q = 1/2*a.^2.*b.^2 - 1/2*d.^2.*a.^2 - 3/2.*b.^2.*c.^2 + 2*a.*b.*c.*d - 1/2.*b.^4.*c.^2./d.^2;
%                 s = d.^2.*c.*a - d.*c.^2.*b + d.*a.^2.*b - a.*b.^2.*c;
%                 fac1 = (d.*SJ + c) ./d;
%                 
%                 
%                 f1 = b.^2.*SJ.^2 + d.^2.*SJ.^2 + 2*a.*b.*SJ + 2*c.*d.*SJ + a.^2 + c.^2;
%                 g1 = (SJ + a.*b + c.*d)./t;
%                 
%         
%                 A1 =  -(xpmat(i,j) - xnew(1:end-1)).*cos(thetpan) - (ypmat(i,j) - ynew(1:end-1)).*sin(thetpan);
%         B1 =  (xpmat(i,j) - xnew(1:end-1)).^2 + (ypmat(i,j) - ynew(1:end-1)).^2 ;
%         E1 =  (xpmat(i,j) - xnew(1:end-1)).*sin(thetpan) - (ypmat(i,j) - ynew(1:end-1)).*cos(thetpan);
%        
%         temp1 = 2*A1.*SJ + SJ.^2 + B1;
%         
%                  ga = atan2(E1.*SJ,(B1+SJ.*A1)); 
%                 gb = log(temp1./B1);
%                 gc =  atan2(a + b.*SJ, c + d.*SJ);
%           % this condition is different when doing on the surface -- bottom of surface has one branch, top another     
%          % oneloc = min(find(xnew<=xpmat(i,j)));
%            % twoloc = max(find(xnew<xpmat(i,j)));   
%              %  if (twoloc < length(ynew))
%                  %   if (ypmat(i,j) <  ( ynew(twoloc) + ynew(twoloc+1))/2)
%              
%                         
%                   yonchord = (ynew(locxmax) - ynew(locxmin))/(xnew(locxmax)-xnew(locxmin))*(xpmat(i,j) - xnew(locxmin)) + ynew(locxmin);
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
%                  else 
%                          for kk = 1:length(ynew)-1
%                     if (a(kk) + b(kk).*SJ(kk) < 0) 
%                         if ( c(kk) + d(kk).*SJ(kk)  < 0 )
%                             gc(kk) = gc(kk) + 2*pi;
%                         end
%                     end 
%                          end
%                     end
%                     
% %                 else         
% %                 
% %                      for kk = 1:length(ynew)-1
% %                     if (a(kk) + b(kk).*SJ(kk) < 0) 
% %                         if ( c(kk) + d(kk).*SJ(kk)  < 0 )
% %                             gc(kk) = gc(kk) + 2*pi;
% %                         end
% %                     end 
% %                     end
% %                     
% %                 end
% %             
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
%            %    for jj = 1:length(fac2)
%               
%                     
%              %   if     (  abs(fac2(jj)) < 1e-8)
% %                 if (abs(t(jj)) < 1e-12)
% %                     
% %                     f3(jj) = 0;
% %                     f4(jj) = 0;
% %                    
% %                 else 
% %                    
% %                     f3(jj) = fac2(jj).*atan(g1(jj))./t(jj);
% %                      f4(jj) = 2*atan(g1(jj)).*fac3(jj)./t(jj);
% %                 end
% %                end
%                
% 
% for jj = 1:length(ztop)
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
%                 
%                 psiint1_SJ = log(f1).*(SJ + a.*b + c.*d) + f4 -2*SJ  ;
%                 psiint2_SJ = 1/2*log(f1).*(SJ.^2 - a.^2.*b.^2+a.^2.*d.^2-4*a.*b.*c.*d+b.^2.*c.^2-c.^2.*d.^2) + f3 ...
%                      - 1/2*SJ.^2 + a.*b.*SJ + c.*d.*SJ;
%                 
%           
%                
%                  
%                 f1 =  a.^2 + c.^2;
%                 g1 = (a.*b + c.*d)./t;
%                 
%                
%                  
%                for jj = 1:length(fac2)
%                 
%                 if (abs(t(jj)) < 1e-12)
%                     
%                     f3(jj) = 0;
%                     f4(jj) = 0;
%                    
%                 else 
%                    
%                     f3(jj) = fac2(jj)*atan(g1(jj))./t(jj);
%                      f4(jj) = 2*atan(g1(jj)).*fac3(jj)./t(jj);
%                 
%                end
%                end
%                
%            
%                 
%                  psiint1_0 = log(f1).*(a.*b + c.*d) + f4  ;
%                 psiint2_0 = 1/2*log(f1).*(- a.^2.*b.^2+a.^2.*d.^2-4*a.*b.*c.*d+b.^2.*c.^2-c.^2.*d.^2) + f3 ;
%                      
%                 
%                
%                 psiansint1 = psiint1_SJ - psiint1_0;
%                 psiansint2 = psiint2_SJ - psiint2_0;
%                        
%                 
%                 anspsipan = gamans(1:end-1)'.*(psiansint1)/4/pi  ...
%                    +    ( gamans(2:end)' - gamans(1:end-1)'  )./SJ.*(psiansint2)/4/pi  ;
%                 
%                 
%             
%                 
%                 tatanup = atan2(ypmat(i,j) ,(xpmat(i,j)+1)  ) ;
%                 tatandown = atan2(ypmat(i,j),(xpmat(i,j)-1) ) ;
%                 
%                 if (ypmat(i,j) > 0)
%                         if (xpmat(i,j)-1 > 0)  
%                             tatandown = tatandown-pi;
%                         end
%                         if (xpmat(i,j)-1 < 0  ) 
%                         tatandown = tatandown - pi; 
%                         end
%                 elseif (ypmat(i,j) < 0) 
%                         if (xpmat(i,j)-1 > 0)
%                             tatandown = tatandown +pi;
%                         end
%                        if (xpmat(i,j) -1 < 0) 
%                             tatandown = tatandown + pi;
%                         end
%                 end
%                     
%                
%       %%%%
%       % new phi way - totally integrated
%       %%%%
%                  phiI = gc.*SJ  + E1/2.*gb   -A1.* ga;
%                 
%                 
%                
%                phiII_overS = E1/2 -ga.*(-2*A1.^2 + B1)./SJ/2 - E1.*A1.*gb./SJ/2 + SJ/2.*gc;
%                
%                 
%                 ansphipan = gamans(1:end-1)'.*(phiI)/2/pi  ...
%                    +    ( gamans(2:end)' - gamans(1:end-1)'  ).*(phiII_overS)/2/pi  ;
%                 philocpannew(i,j) = - sum(real(ansphipan));
%                 philocnew(i,j) = -sum(real(ansphipan))   + delta/4/pi*(  log( (xpmat(i,j)+1)^2+ypmat(i,j)^2) - log( (xpmat(i,j)-1)^2 + ypmat(i,j)^2) ) ...
%                       +   gup/2/pi* tatanup  + gamdown/2/pi*tatandown ;
%         %%%%%%%%%%%%%                
%                 
%                 psilocpan(i,j) = sum(anspsipan);
%                 psiloc(i,j) = sum(anspsipan) +delta/2/pi*(tatanup - tatandown) ...
%                        - gup/4/pi*log( (xpmat(i,j) +1)^2 + ypmat(i,j)^2)  -gamdown/4/pi*log( (xpmat(i,j)-1)^2 + ypmat(i,j)^2);
%                     
%                 
%             end
%         end
%        z2 = xpmat + 1i*ypmat;
%         z1 = atanh(z2)*h/pi;
%         %z1 = xpinaf + 1i*ypinaf;
%         dz2dz1 = pi/h*(1 - (tanh(pi*z1/h)).^2 );
%         compv2 = upan - 1i*vpan;
%         compv1 = compv2.*dz2dz1; 
%             figure(22)
%         plot(real(z1),psiloc)
%         hold on
%         plot(real(z1),philocnew,'r');
        
        