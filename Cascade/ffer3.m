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
clear f fer2
vortstrength=-0.02;
izeropad =1;


%D = BEM2d_readbin('2example');

% time = 0.0125:0.0125:10;
% time = 0.01:0.01:8.0;
%foft = lift_coef;

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
fer = fft(foft,izeropad*length(foft));

hlen = length(fer)/2;
for i=1:hlen,
fer2(i) = fer(i)/numhalf;
f(i)=(i-1)/numpoints*fs;
end;
figure(1);
hh1 = subplot(3,1,2);
%set(hh1,'fontsize',14);
semilogy(f,abs(fer2)/vortstrength);
axis([0 7 -.3 .3]);
hh2 = subplot(3,1,3);
plot(f,imag(fer2)/vortstrength);
axis([0 7 -.3 .3]);
