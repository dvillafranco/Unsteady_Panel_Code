%  Program takes the output of the BEM code 
%	1) plots it 
%	2) transforms it into the frequency domain
%	3) Plots its transform
%
%  user has to put in vortstrength or have the
%     program read it from Vort.in
%  user should make sure that the lift has a portion of 
%     in which the steady state lift is achieved and the
%     vortex has not yet reached the airfoil.
%	one can use wakemovie.m for this.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

vortstrength=0.01;

%D = BEM2d_readbin;
D = 1;
time = 0.0:0.001:0.168;

%%%
% try to position the peak in lift (when vortex influence is
%	first felt)  at the beginning of the time series of data
%   right now I am finding the peak and moving back from it 
%	3 times steps.  This number of times steps may need to 
%	be based on the dt or a similar parameter.
%%%%
% [xmax,iloc]=max(D.L);
% numsteady=1;
% time = D.t(numsteady:length(D.t));
% %foft = D.L(numsteady:length(D.t))-D.L(numsteady);
 foft = lift_coef*10^14;
  %foft = test_var;

dt = time(2) - time(1);
fs = 1/dt;
fshalf = fs/2;
%
% number of data points in file and half it
%	is nyquest number past that is 
%	mirror of negative frequencies
%
numpoints = length(time);
numhalf = numpoints/2;
%
% data file with the data to be transformed
%
f1 =figure(1);
set(f1,'Position',[25 180 412 484]);
hh1 = subplot(3,1,1);
plot(time,foft);
%
% transform it, only keep the first half and
%	get ready to plot on real time and
%	amplitude scale
%
fer = fft(foft);

hlen = length(fer)/2;
for i=1:hlen,
fer2(i) = fer(i)/numhalf;
f(i)=(i-1)/numpoints*fs;
end;
figure(1);
hh1 = subplot(3,1,2);
%set(hh1,'fontsize',14);
plot(f,(real(fer2).^2+imag(fer2).^2).^0.5/vortstrength);
%axis([0 f(hlen) -.1 .3]);
%axis([0 30 -.1 .3]);

hh2 = subplot(3,1,3);
plot(f,imag(fer2)/vortstrength);
%axis([0 f(hlen) -.1 .3]);
%axis([0 30 -.1 .3]);

figure (2)
plot(f,angle(fer2)/vortstrength);

figure(3)
plot(f,(real(fer2).^2+imag(fer2).^2).^0.5,'LineWidth',2)
set(gca,'FontSize',16)
ylabel('Transform of Lift')
xlabel('F')

