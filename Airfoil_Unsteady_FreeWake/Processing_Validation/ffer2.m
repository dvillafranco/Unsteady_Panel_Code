%  Program takes a sine wave in time 
%	1) plots it 
%	2) transforms it into the frequency domain using ZERO PADDING
%	3) Plots its transform
%  Does the inverse transform assuming that only the
%	real data (1/2 of the FFT) is known
%	4) Sets up frequency domain series
%	5) Plots it
%	6) Does the IFFT
%	7) Plots it
%  One can just feed the output of the FFT back into an
% 	IFFT and get the right solution!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

vortstrength=0.02;
izeropad =3.0;
% time = -4.50:0.01:4.49;
time = 0.01:0.01:30;

%time =parttime; 
%foft = partlift;
%time = D.t; 
%foft= [zeros(1,length(lift_coef)),lift_coef];
foft = lift_coef;

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
fer2(i) = fer(i)/(numhalf);
f(i)=(i-1)/length(fer)*fs;
end;
figure(1);
hh1 = subplot(3,1,2);
%set(hh1,'fontsize',14);
plot(f,real(fer2)/vortstrength);
plot(f,abs(fer2)/vortstrength);
%axis([0 20 -.3 .3]);
hh2 = subplot(3,1,3);
%plot(f,imag(fer2)/vortstrength);
plot(f,unwrap(angle(fer2))/vortstrength);
%axis([0 20 -.3 .3]);

%--- supposedly method for obtaining power spectral density
% listed in FFT in matlab help
%Pyy = fer.*conj(fer) /length(time);
%plot(f, Pyy(1:length(fer2)));



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% now try to do the inverse fft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

df = f(2) - f(1);
ts = 1/df;
tshalf = ts/2;
numpoints2 = length(fer2);
numhalf2 = numpoints2/2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% have to recover time frequency domain signal as if 
%	it came from an fft in matlab
% If you are only performing an ifft from freq. data instead of
% 	getting the freq. data from an fft in matlab you
%	will not know numhalf.
% Thus you just have to use length(fer2) and hope that the
%	original function did not come from something with
%	a zero padding.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
signal = fer2.*numhalf;
%
f2=figure(3);
set(f2,'Position',[525 230 412 414]);
hh1 = subplot(2,1,1);
plot(f,real(signal),f,imag(signal),'--');
%
% 
% This function generates the inverse Fourier Transform
% from only the positive frequency complex FFT vector.
% E.J. O'Neill
% 5 Nov. 1996
% Background:
% The Matlab convention that must be understood to generate
% this program is that for an even number N point FFT,
% point 1 is the DC component and point (N./2)+1 is zero.
% Therefore, from only the positive frequency component, the
% time series can be generated by performing the IFFT on a
% vector which is 
%full_fft=[posfft(1....N) 0 posfft(N....2)]
%[m,n]=size(posf_fft);

% if data comes from fft from mat lab this works
% otherwise need neg_freq=[signal(2:numpoints2).' 0];
% and pos_freq = signal.'
%
neg_freq=[signal(2:numpoints2) 0];
neg_freq=fliplr(neg_freq);
pos_freq=signal;
neg_freq=conj(neg_freq);
fullfft=[pos_freq neg_freq];
%ts=real(ifft(fullfft).');

siner1 = (ifft(fullfft).');

time1 = (0:length(siner1)-1)*ts/length(siner1);
%end;
figure(3);
hh1 = subplot(2,1,2);
set(hh1,'fontsize',14);
plot(time1,real(siner1));
%axis([0,ts,-1,1]); 
%axis([0,time(end)+2,-.1,.1]); 
