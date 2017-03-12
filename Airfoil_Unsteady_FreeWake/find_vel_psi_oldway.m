% get velocity field around object  (currently object goes from -1 to 1 and
% then  -.2  to .2 

 % xp =  [-3.801 : .2 :4.8 ];
  xp =  [-.61 : .02 : .6 ];
  yp = [ -.21 : .02 : .2];
 % yp = [-1.5 : .1 : 1.5  ]%    .025:.01:.2];
 
% above airfoil
% xp = [ -.701 : .01 : .7];
 %yp = [ -.06  -.05  .05 .06];
  
  
  
numx = length(xp);
numy = length(yp);
xpmat = [];
for i = 1:numy
       xpmat = [xpmat ; xp];
end
ypmat = [];
for j = 1:numx
    ypmat = [ypmat yp'];
end



xpmat = xmidmat;
ypmat = ymidmat;


numx = length(xmid);
numy = length(ymid);



for j = 1:numx
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
           
           up = gamans(1:end-1)'/2/pi.* (  (-ypmat(i,j)+ynew(1:end-1)).*(ga) + sin(thetpan).*(gb)  )+ ...
               (gamans(2:end)'-gamans(1:end-1)')./2/pi./SJ.* ( ( -ypmat(i,j) + ynew(1:end-1)).*(gb) + sin(thetpan).*(gc));
           upan(i,j) = -sum(up);
           
           upj = gamans(1:end-1)'/2/pi.* (  (-ypmat(i,j)+ynew(1:end-1)).*(ga) + sin(thetpan).*(gb)  )+ ...
               (-gamans(1:end-1)')./2/pi./SJ.* ( -ypmat(i,j) + ynew(1:end-1)).*(gb) + sin(thetpan).*(gc);
           upjp1 = (gamans(2:end)')./2/pi./SJ.* ( -ypmat(i,j) + ynew(1:end-1)).*(gb) + sin(thetpan).*(gc);
           
           
           vp = gamans(1:end-1)'/2/pi.* (  (xpmat(i,j) - xnew(1:end-1)).*(ga) - cos(thetpan).*(gb) ) + ...
               (gamans(2:end)'-gamans(1:end-1)') ./2/pi./SJ.*  (  (xpmat(i,j) - xnew(1:end-1)).*(gb) - cos(thetpan).*(gc));
          
           
           vpj =  gamans(1:end-1)'/2/pi.* (  (xpmat(i,j) - xnew(1:end-1)).*(ga) - cos(thetpan).*(gb) ) + ...
               (-gamans(1:end-1)') ./2/pi./SJ.*  (  (xpmat(i,j) - xnew(1:end-1)).*(gb) - cos(thetpan).*(gc));
           vpjp1 = (gamans(2:end)') ./2/pi./SJ.*  (  (xpmat(i,j) - xnew(1:end-1)).*(gb) - cos(thetpan).*(gc));
          
           
           vpan(i,j) = -sum(vp);
           
%          upj*(-sin(thetpan(i))) + vpj*(cos(thetpan(i)));
%          gamj(i,j) ;
           
     
          
           % get source, sink and vortices contribution to velocity field
           
           upan(i,j) = upan(i,j) + vinf*cos(alp);
           
           vpan(i,j) = vpan(i,j) + vinf * sin(alp);
          
    end
end
    figure(1)
        hold on;
    
        quiver( xpmat, ypmat, upan, vpan,1)  
        
    cp = 1 -(upan.^2 + vpan.^2); 
    figure(5)
    plot(xpmat,cp);

    % when computing velocity field on body - use this to check normal
    % condition
    
%         for i = 1:length(xmid)
%             checker = upan(i,i)*(-sin(thetpan(i))) + vpan(i,i)*cos(thetpan(i)) 
%         end;
%         
[maxxnew locxmax] = max(xnew);
[minxnew locxmin] = min(xnew);
        
        %Build potential field due to panels
        
        % Loop over points where you want the value of the potential
        for j = 1:numx
            for i = 1:numy
             
                
                a = ypmat(i,j) - ynew(1:end-1);
                b = -sin(thetpan);
                c = xpmat(i,j) - xnew(1:end-1);
                d = -cos(thetpan);
                t = (a.*d-b.*c);
                
                
                f1 = b.^2.*SJ.^2 + d.^2.*SJ.^2 + 2*a.*b.*SJ + 2*c.*d.*SJ + a.^2 + c.^2;
                g1 = (SJ + a.*b + c.*d)./t;
                
  % repeat terms as in the velocity part to help with the potential
  % function
  
   A1 =  -(xpmat(i,j) - xnew(1:end-1)).*cos(thetpan) - (ypmat(i,j) - ynew(1:end-1)).*sin(thetpan);
        B1 =  (xpmat(i,j) - xnew(1:end-1)).^2 + (ypmat(i,j) - ynew(1:end-1)).^2 ;
        E1 =  (xpmat(i,j) - xnew(1:end-1)).*sin(thetpan) - (ypmat(i,j) - ynew(1:end-1)).*cos(thetpan);
       
        temp1 = 2*A1.*SJ + SJ.^2 + B1;
        
                 ga = atan2(E1.*SJ,(B1+SJ.*A1)); 
                gb = log(temp1./B1);
                gc =  atan2(a + b.*SJ, c + d.*SJ);
              
                
           
                
%                 %  y - ymin = (ymax - ymin)/(xmax-xmin)*(x - xmin)
                        yonchord = (ynew(locxmax) - ynew(locxmin))/(xnew(locxmax)-xnew(locxmin))*(xpmat(i,j) - xnew(locxmin)) + ynew(locxmin);
                       
                if ( xpmat(i,j) > minxnew & xpmat(i,j) < maxxnew & ypmat(i,j) < yonchord) 
                    
                    
                    for kk = 1:length(ynew)-1
                    if (a(kk) + b(kk).*SJ(kk) < 0) 
                            gc(kk) = gc(kk) + 2*pi;
                    end 
                    if (a(kk) + b(kk).*SJ(kk) > 0 & c(kk) + d(kk).*SJ(kk) > 0) 
                        gc(kk) = gc(kk) + 2*pi;
                    end
                    end
                    
                else   
                    
                for kk = 1:length(ynew)-1
                    if (a(kk) + b(kk).*SJ(kk) < 0) 
                        if ( c(kk) + d(kk).*SJ(kk)  < 0 )
                            gc(kk) = gc(kk) + 2*pi;
                        end
                       
                    end
                    
                end
                end
                                
                
               
                
                fac2 = (-a.^3.*b-a.^2.*c.*d-a.*b.*c.^2-c.^3.*d+a.^3.*b.^3-a.^3.*b.*d.^2+5*a.^2.*b.^2.*c.*d - a.^2.*c.*d.^3 - a.*b.^3.*c.^2 ...
                    + 5.*a.*b.*c.^2.*d.^2 - b.^2.*c.^3.*d + c.^3.*d.^3);
      
                fac3 = (a.^2 + c.^2 - a.^2.*b.^2 - 2*a.*b.*c.*d - c.^2.*d.^2);
               
                for jj = 1:length(fac2)
                
                    
             %   if     (  abs(fac2(jj)) < 1e-8)
                if (abs(t(jj)) < 1e-12)
                    
                    f3(jj) = 0;
                    f4(jj) = 0;
                   
                else 
                   
                    f3(jj) = fac2(jj).*atan(g1(jj))./t(jj);
                     f4(jj) = 2*atan(g1(jj)).*fac3(jj)./t(jj);
                end
               end
                
                
                psiint1_SJ = log(f1).*(SJ + a.*b + c.*d) + f4 -2*SJ  ;
                psiint2_SJ = 1/2*log(f1).*(SJ.^2 - a.^2.*b.^2+a.^2.*d.^2-4*a.*b.*c.*d+b.^2.*c.^2-c.^2.*d.^2) + f3 ...
                     - 1/2*SJ.^2 + a.*b.*SJ + c.*d.*SJ;
                
          
               
                 
                f1 =  a.^2 + c.^2;
                g1 = (a.*b + c.*d)./t;
                
                
                
      
                 
               for jj = 1:length(fac2)
                  
                   
              
                if (abs(t(jj)) < 1e-12)
                    
                    f3(jj) = 0;
                    f4(jj) = 0;
                   
                else 
                   
                    f3(jj) = fac2(jj)*atan(g1(jj))./t(jj);
                     f4(jj) = 2*atan(g1(jj)).*fac3(jj)./t(jj);
                
               end
               end
        
                 psiint1_0 = log(f1).*(a.*b + c.*d) + f4  ;
                psiint2_0 = 1/2*log(f1).*(- a.^2.*b.^2+a.^2.*d.^2-4*a.*b.*c.*d+b.^2.*c.^2-c.^2.*d.^2) + f3 ;
                     
       
                
                psiansint1 = psiint1_SJ - psiint1_0;
                psiansint2 = psiint2_SJ - psiint2_0;
                
      %%%%
      % new phi way - totally integrated
      %%%%
                 phiI = gc.*SJ  + E1/2.*gb   -A1.* ga;
                
                
               
               phiII_overS = E1/2 -ga.*(-2*A1.^2 + B1)./SJ/2 - E1.*A1.*gb./SJ/2 + SJ/2.*gc;
               
                
                ansphipan = gamans(1:end-1)'.*(phiI)/2/pi  ...
                   +    ( gamans(2:end)' - gamans(1:end-1)'  ).*(phiII_overS)/2/pi  ;
                philocpannew(i,j) = - sum(real(ansphipan));
                philocnew(i,j) = -sum(real(ansphipan))  + vinf*cos(alp)*xpmat(i,j) + vinf*sin(alp)*ypmat(i,j);
        %%%%%%%%%%%%%    
              
                anspsipan = gamans(1:end-1)'.*(psiansint1)/4/pi  ...
                   +    ( gamans(2:end)' - gamans(1:end-1)'  )./SJ.*(psiansint2)/4/pi  ;
                
                
                psilocpan(i,j) = sum(anspsipan);
                psiloc(i,j) = sum(anspsipan) + vinf*cos(alp)*ypmat(i,j) - vinf*sin(alp)*xpmat(i,j);
                

            end
        end
        
        figure(10)
        surf(xpmat,ypmat,(philocnew))
        
         figure(12)
        surf(xpmat,ypmat,(psiloc))
        
        figure(13)
        pcolor(xpmat,ypmat,psiloc);
        hold on
        plot(xnew,ynew);
        shading interp;
        
        cs = [-.25 :.02 :.25];
        cs2 = [-2 :.08 :2];
        figure(14);
        contour(xpmat,ypmat,psiloc,40);
        hold on
     %   contour(xpmat,ypmat,philocnew,40);
        plot(xnew,ynew)
        
        % check potential function    dphi / dx = u     dphi / dy = v
        
        dphidxbot = (philocnew(2,2:end) - philocnew(2,1:end-1))/.05;
        dphidxtop = (philocnew(end-1,2:end) - philocnew(end-1,1:end-1))/.05;
        dphidybot = (philocnew(2,:) - philocnew(1,:))/.01;
        dphidytop = (philocnew(end-2,:) - philocnew(end-3,:))/.01;
        
        figure(20);
        plot(xp(1:end-1), dphidxbot, xp(1:end), upan(2,:));
        title('u bottom');
        figure(21);
        plot(xp(1:end-1), dphidxtop, xp(1:end), upan(end-1,:));
            title('u top');
        figure(22);
        plot(xp(1,:), dphidybot, xp(1:end), vpan(2,:));
        title('v bottom');
        figure(23);
        plot(xp(1,:), dphidytop, xp(1:end), vpan(end-2,:)); 
        title('v top');