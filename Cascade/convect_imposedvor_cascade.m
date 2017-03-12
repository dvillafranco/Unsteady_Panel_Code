clear xpmat ypmat upan vpan psiloc philoc  psilocpan philocpan philocnew

 
xp = x_imp(n);
yp = y_imp(n);

numx = length(xp);
numy = length(yp);

for i = 1:numy
    xpmat = xp;
end

for j = 1:numx
    ypmat = yp';
end


bot33 = (xpmat - 1).^2 + ypmat.^2;
bot44 = (xpmat + 1).^2 + ypmat.^2;

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
           
          
           
           %gamdown = sum( gamans(2:end-1)'.*SJ/2 + gamans(1:end-2)'.*SJ/2) - vinf*h*sin(alp);
           
           % get source, sink and vortices contribution to velocity field
           
           upan_imp(i,j) = upan(i,j) + delta/2/pi*  ((xpmat(i,j)+1)./bot44(i,j) - (xpmat(i,j)-1)./bot33(i,j) ) -  gup /2/pi*(ypmat(i,j)./bot44(i,j)) - gamdown/2/pi*(ypmat(i,j)./bot33(i,j))+...
               sum((gam_dim./(1)).*((ypmat(i,j)-y_fixed(1:end-1))./((xpmat(i,j)-x_fixed(1:end-1)).^2+(ypmat(i,j)-y_fixed(1:end-1)).^2)));
           
           vpan_imp(i,j) = vpan(i,j) + delta/2/pi*(  ypmat(i,j)./bot44(i,j) - ypmat(i,j)./bot33(i,j)) + gup/2/pi*( (xpmat(i,j)+1)./bot44(i,j) ) + gamdown/2/pi*((xpmat(i,j)-1)./bot33(i,j))+ ...
               sum((-gam_dim/(1)).*((xpmat(i,j)-x_fixed(1:end-1))./((xpmat(i,j)-x_fixed(1:end-1)).^2+(ypmat(i,j)-y_fixed(1:end-1)).^2)));
           
    end
end


%dt = 0.01;
x_imp(nx+1) = x_imp(nx)+upan_imp*dt;
y_imp(nx+1) = y_imp(nx)+vpan_imp*dt;


% z2 = xpmat + 1i*ypmat;
%         z1 = atanh(z2)*h/pi;
%         %z1 = xpinaf + 1i*ypinaf;
%         dz2dz1 = pi/h*(1 - (tanh(pi*z1/h)).^2 );
%         compv2 = upan - 1i*vpan;
%         compv1 = compv2.*dz2dz1;
%         
%         u1 = real(compv1);
%         v1 = -imag(compv1);
%             
%        figure(2)
%         hold on;
%        % quiver(xpinaf,ypinaf, u1, v1);
%         
%         quiver(real(z1),imag(z1), u1, v1);
%         
% %         
% %         for i = 1:length(xmid)
% %             checker = upan(i,i)*(-sin(thetpan(i))) + vpan(i,i)*cos(thetpan(i)) 
% %         end;
% 
%         figure(5)
%         quiver( xpmat, ypmat, upan, vpan,6)
%         
%         figure(6)
%         quiver(real(z1),imag(z1), u1, v1,3);
%         
%   
% lastpoints =  round(length(ynew)/2); 
%  [maxxnew locxmax] = max(xnew);
% [minxnew locxmin] = min(xnew);
