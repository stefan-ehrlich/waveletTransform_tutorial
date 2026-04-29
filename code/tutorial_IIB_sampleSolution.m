%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Tutorial II part B - sample solution
% Author: Stefan Ehrlich
% Date: 19.01.2020
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

inpath = strcat(pwd,'\data\');
load(strcat(inpath,'dummyX.mat'));


%% 1. Convolution in the time- and frequency domain
signal = dummyX.x;
kernel = dummyX.gausskernel;

halfkernel = (length(kernel)-1)/2;


%% 1a. compute convolution via MATLAB 'conv' function
conv_result1 = conv(signal,kernel); % convolve signal and kernel
conv_result1 = conv_result1(halfkernel+1:end-halfkernel); % prune convolved data

%% 1b. compute convolution via piecewise dot product
signal_pad = [zeros(1,halfkernel) signal zeros(1,halfkernel)]; % zero padding signal

for w = 1:length(signal_pad)-length(kernel)+1
    conv_result2(w) = sum(signal_pad(w:w+length(kernel)-1).*kernel);
end

%% 1c. compute convolution via multiplication of spectra in frequency domain
nConv = length(signal)+length(kernel)-1; % N+K-1
signal_pad = [signal zeros(1,nConv-length(signal))]; % zero pad signal
kernel_pad = [kernel zeros(1,nConv-length(kernel))]; % zero pad kernel

SIGNAL = fft(signal_pad); % compute FFT of zero-padded signal
KERNEL = fft(kernel_pad); % compute FFT of zero-padded kernel

RESULT3 = SIGNAL.*KERNEL; % multiply in frequency domain

conv_result3 = ifft(RESULT3); % inverse FFT
conv_result3 = conv_result3(halfkernel+1:end-halfkernel); % prune resulting data


figure(1)
subplot(3,1,1)
plot(signal)
title('Signal')
grid on
subplot(3,1,2)
plot(-20:1:20,kernel)
title('Gauss-Kernel')
grid on
subplot(3,1,3)
plot(conv_result1)
hold on
plot(conv_result2)
plot(conv_result3)
grid on
legend('MATLAB conv','manual conv.','fft-based conv.')

%% 2. Create a single wavelet
% kernel-length = 1 sec; f = 10 Hz; nCycles = 5; % length of kernel = 1 sec
kernel_length = 1; 
wavtime = -kernel_length/2:1/dummyX.srate:kernel_length/2;

nCycles = 5;
f=10; s=nCycles/(2*pi*f); %A=1/((s*sqrt(pi)).^(1/2));
gaussWin = exp((-wavtime.^2)/(2*s^2));
%gaussWin = gaussWin./sum(gaussWin);
complSine = exp(1i*2*pi*f.*wavtime);
cmw = gaussWin.*complSine;
halfwave = floor(length(cmw)/2)+1;

figure(2)
subplot(1,2,1)
plot(wavtime,gaussWin)
hold on
plot(wavtime,real(complSine))
plot(wavtime,real(cmw),'k','LineWidth',1.5)
xlabel('time [ms]')
subplot(1,2,2)
plot3(wavtime,real(cmw),imag(cmw),'k','LineWidth',1.5)
grid on
title('Complex Morlet Wavelet')
xlabel('time [ms]')
ylabel('real')
zlabel('imag')

%% 3. Complex wavelet transform with dummy signal
signal = dummyX2.x;
% zero padding
nConv = length(signal)+length(cmw)-1;
signal_pad = [signal zeros(1,nConv-length(signal))];
cmw_pad = [cmw zeros(1,nConv-length(cmw))];

% Frequency domain convolution
SIGNAL = fft(signal_pad);
CMW = fft(cmw_pad);
CMW = CMW./max(CMW); % normalize kernel FFT
RESULT = SIGNAL.*CMW; % multiply FFT of signal and kernel
cmw_result = ifft(RESULT); % inverse FFT
cmw_result = cmw_result(halfwave-1:end-halfwave);

figure(3)
subplot(3,1,1)
plot(dummyX2.t,dummyX2.x);
title('original signal')
ylabel('Amp.')
grid on
subplot(3,1,2)
plot(dummyX2.t,2*real(cmw_result),'LineWidth',1)
hold on
plot(dummyX2.t,2*abs(cmw_result),'LineWidth',1)
plot(dummyX2.t,(2*abs(cmw_result)).^2,'LineWidth',1)
ylabel('Amp. / Amp.²')
legend('bandfiltered signal','magnitude','power')
title('Wavelet transformed data')
grid on
subplot(3,1,3)
plot(dummyX2.t,angle(cmw_result),'LineWidth',1)
xlabel('Time [s]')
ylabel('rad')
title('Phase')
grid on

%% 4. Create and apply complex Morlet wavelet family
% linear spaced frequencies from 4 Hz to 250 Hz in steps of 2 Hz
frex = 4:2:250;
% nCycles from 2 to 10
nCycles = linspace(3,10,length(frex));
% kernel length = 1 sec
kernel_length = 1; 
wavtime = -kernel_length/2:1/dummyX2.srate:kernel_length/2;

% create complex Morlet wavelets

for f=1:length(frex)
    s=nCycles(f)/(2*pi*frex(f)); %A=1/((s*sqrt(pi)).^(1/2));
    cmw_family(f,:) = exp((-wavtime.^2)/(2*s^2)).*exp(1i*2*pi*frex(f).*wavtime);
end

figure(4) 
subplot(1,2,1)
plot(wavtime,real(cmw_family))
grid on
subplot(1,2,2)
imagesc(flip(real(cmw_family)))

for f=1:length(frex)
    cmw_pad = [cmw_family(f,:) zeros(1,nConv-length(cmw_family(f,:)))]; % zero-padding wavelet
    CMW = fft(cmw_pad);
    CMW = CMW./max(CMW); % normalize kernel FFT
    tf = ifft(SIGNAL.*CMW);
    TF(f,:) = tf(halfwave-1:end-halfwave);
    f
end

% visualize
figure(5)
contourf(dummyX2.t,frex,2*abs(TF),10,'linecolor','none')
%set(gca, 'YScale', 'log')
title('CWT Magnitude')
xlabel('Time [s]')
ylabel('frequency [Hz]')
colorbar


