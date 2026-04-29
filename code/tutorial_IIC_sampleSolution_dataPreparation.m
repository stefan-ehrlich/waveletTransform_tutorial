%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Tutorial II part C - sample solution
% Author: Stefan Ehrlich
% Date: 21.01.2020
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

inpath = strcat(pwd,'\data\');
load(strcat(inpath,'DATA_256Hz.mat'));

%% Compute time-frequency data
% linear spaced frequencies from 4 Hz to 250 Hz in steps of 2 Hz
frex = 2:2:40;
% nCycles from 2 to 10
nCycles = linspace(3,10,length(frex));
% kernel length = 1 sec
kernel_length = 2; 
wavtime = -kernel_length/2:1/DATASET.srate:kernel_length/2;
halfwave = floor(length(wavtime)/2)+1;

nConv = size(DATASET.data,2)+length(wavtime)-1;

% create complex Morlet wavelets

for f=1:length(frex)
    s=nCycles(f)/(2*pi*frex(f)); %A=1/((s*sqrt(pi)).^(1/2));
    cmw_family(f,:) = exp((-wavtime.^2)/(2*s^2)).*exp(1i*2*pi*frex(f).*wavtime);
end



for ch = 1:size(DATASET.data,1)
    for tr = 1:size(DATASET.data,3)
        for s = 1:size(DATASET.data,4)
            data = DATASET.data(ch,:,tr,s);
            data_pad = [data zeros(1,nConv-length(data))]; % zero-padding data
            SIGNAL = fft(data,nConv);

            for f=1:length(frex)
                CMW = fft(cmw_family(f,:),nConv);
                CMW = CMW./max(CMW); % normalize kernel FFT
                tf = ifft(SIGNAL.*CMW);
                TF(ch,f,:,tr,s) = tf(halfwave-1:end-halfwave);
            end
        end
    end
    ch
end

TF_DATA = DATASET;
TF_DATA.data = 2*abs(TF);
TF_DATA.chanlocs = DATASET.chanlocs;
TF_DATA.data_dimensions = {'channels','frequencies','time','trials','subjects'};
TF_DATA.frex = frex;

% magnitude computation and baseline correction
for s=1:17
    for tr=1:150
        TF_DATA.data_corr(:,:,:,tr,s)=TF_DATA.data(:,:,:,tr,s)-mean(TF_DATA.data(:,:,TF_DATA.times>-300&TF_DATA.times<-100,tr,s),3);
    end
    s
end
% compute condition average and difference
for s=1:17
    TF_DATA.nonerror(:,:,:,s) = mean(TF_DATA.data_corr(:,:,:,TF_DATA.labels(:,s)==-1,s),4);
    TF_DATA.error(:,:,:,s) = mean(TF_DATA.data_corr(:,:,:,TF_DATA.labels(:,s)==1,s),4);
    TF_DATA.diff(:,:,:,s) = TF_DATA.error(:,:,:,s)-TF_DATA.nonerror(:,:,:,s);
    s
end

% try downsampling data by factor 4
ds = 2:4:512;


tfviewerx(TF_DATA.times(ds),TF_DATA.frex,double(mean(TF_DATA.diff(:,:,ds,:),4)),TF_DATA.chanlocs)

TF_DATA.magnitude = TF_DATA.data_corr(:,:,ds,:,:);
TF_DATA.times = TF_DATA.times(ds);

