%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Tutorial II part A - sample solution
% Author: Stefan Ehrlich
% Date: 12.01.2020
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

load DATA_256Hz;

%% 1. Familiarize with dataset
% XXX

%% 2. Baseline correction
time_bl = DATASET.times<0 & DATASET.times>-200;
DATASET.data = DATASET.data-mean(DATASET.data(:,time_bl,:,:),2);

%% 3. Plotting some data
ch = 20; subject = 1; trial1 = 5; trial2 = 14;

figure(1)
plot(DATASET.times,DATASET.data(ch,:,trial1,subject),'LineWidth',1.5);
hold on
plot(DATASET.times,DATASET.data(ch,:,trial2,subject),'LineWidth',1.5);
xlabel('Time [ms]'); ylabel(DATASET.data_unit)
title('Channel Cz, trial 5 and 14')
legend('non-error','error')
grid on

%% 4. Single channel ERP grand average time-course
for s=1:length(DATASET.subjects)
    ERP.nonerror_avg(:,:,s) = mean(DATASET.data(:,:,DATASET.labels(:,s)==-1,s),3);
    ERP.error_avg(:,:,s) = mean(DATASET.data(:,:,DATASET.labels(:,s)==1,s),3);
end
ERP.diff_avg = ERP.error_avg-ERP.nonerror_avg;
ERP.nonerror_Gavg = mean(ERP.nonerror_avg,3)
ERP.error_Gavg = mean(ERP.error_avg,3)
ERP.diff_Gavg = ERP.error_Gavg-ERP.nonerror_Gavg;

sel_time = DATASET.times>-200 & DATASET.times<1200;

figure(2)
plot(DATASET.times(sel_time),ERP.nonerror_Gavg(ch,sel_time),'b','LineWidth',1.5)
hold on
plot(DATASET.times(sel_time),ERP.error_Gavg(ch,sel_time),'r','LineWidth',1.5)
plot(DATASET.times(sel_time),ERP.diff_Gavg(ch,sel_time),'k','LineWidth',1.5)
title('timecourse over channel Cz')
xlabel('time locked to feedback [ms]')
ylabel('Amplitude [µV]')
legend('non-error','error','diff: error minus non-error')
grid on



%% 3xxx. Global Field Power (GFP) and Global Map Dissimilarity (DISS)
% single subject GFP
singleGFP.nonerror_Gavg = mean(mean(ERP.nonerror_avg.^2,1).^(1/2),3);
singleGFP.error_Gavg = mean(mean(ERP.error_avg.^2,1).^(1/2),3);
singleGFP.diff_Gavg = singleGFP.error_Gavg-singleGFP.nonerror_Gavg;

% grand average GFP
grandGFP.nonerror_Gavg = mean(ERP.nonerror_Gavg.^2,1).^(1/2);
grandGFP.error_Gavg = mean(ERP.error_Gavg.^2,1).^(1/2);
grandGFP.diff_Gavg = grandGFP.error_Gavg-grandGFP.nonerror_Gavg;

figure(3)
subplot(2,1,1)
plot(DATASET.times(sel_time),singleGFP.nonerror_Gavg(sel_time),'b','LineWidth',1.5)
hold on
plot(DATASET.times(sel_time),singleGFP.error_Gavg(sel_time),'r','LineWidth',1.5)
plot(DATASET.times(sel_time),singleGFP.diff_Gavg(sel_time),'k','LineWidth',1.5)
title('Global Field Power (GFP) - single subject average')
xlabel('time locked to feedback [ms]')
ylabel('Power [µV²]')
legend('non-error','error','diff: error minus non-error')
grid on
subplot(2,1,2)
plot(DATASET.times(sel_time),grandGFP.nonerror_Gavg(sel_time),'b','LineWidth',1.5)
hold on
plot(DATASET.times(sel_time),grandGFP.error_Gavg(sel_time),'r','LineWidth',1.5)
plot(DATASET.times(sel_time),grandGFP.diff_Gavg(sel_time),'k','LineWidth',1.5)
title('Global Field Power (GFP) - grand average')
xlabel('time locked to feedback [ms]')
ylabel('Power [µV²]')
legend('non-error','error','diff: error minus non-error')
grid on

%% 5. Global Map Dissimilarity
% Normalize maps to unit global field power
% single subject GMD
for s=1:17
    singleGMD.nonerror_avg(:,:,s) = ERP.nonerror_avg(:,:,s)./(mean(ERP.nonerror_avg(:,:,s).^2,1).^(1/2));
    singleGMD.error_avg(:,:,s) = ERP.error_avg(:,:,s)./(mean(ERP.error_avg(:,:,s).^2,1).^(1/2));
end
singleGMD.diff_Gavg = mean(mean((singleGMD.error_avg-singleGMD.nonerror_avg).^2,1).^(1/2),3)

grandGMD.nonerror_avg = ERP.nonerror_Gavg./(mean(ERP.nonerror_Gavg.^2,1).^(1/2));
grandGMD.error_avg = ERP.error_Gavg./(mean(ERP.error_avg.^2,1).^(1/2));

grandGMD.diff_Gavg = mean((grandGMD.error_avg-grandGMD.nonerror_avg).^2,1).^(1/2);

figure(4)
subplot(2,1,1)
plot(DATASET.times(sel_time),singleGMD.diff_Gavg(sel_time),'b','LineWidth',1.5)
title('Global Map Dissimilarity (GMD) - single subject average')
xlabel('time locked to feedback [ms]')
ylabel('unitless')
grid on

subplot(2,1,2)
plot(DATASET.times(sel_time),grandGMD.diff_Gavg(sel_time),'b','LineWidth',1.5)
title('Global Map Dissimilarity (GMD) - grand average')
xlabel('time locked to feedback [ms]')
ylabel('unitless')
grid on

% overview ERP (Cz), GFP difference, GMD
figure(5)
subplot(2,1,1)
plot(DATASET.times(sel_time),zscore(ERP.diff_Gavg(20,sel_time)),'k','LineWidth',1)
hold on
plot(DATASET.times(sel_time),zscore(singleGFP.diff_Gavg(sel_time)),'b','LineWidth',1.5)
plot(DATASET.times(sel_time),zscore(singleGMD.diff_Gavg(sel_time)),'m','LineWidth',1.5)
title('Overview')
xlabel('time locked to feedback [ms]')
ylabel('z-score normalized')
legend('ERP @ Cz','grand GFP','grand GMD')
grid on

subplot(2,1,2)
plot(DATASET.times(sel_time),zscore(ERP.diff_Gavg(20,sel_time)),'k','LineWidth',1)
hold on
plot(DATASET.times(sel_time),zscore(grandGFP.diff_Gavg(sel_time)),'b','LineWidth',1.5)
plot(DATASET.times(sel_time),zscore(grandGMD.diff_Gavg(sel_time)),'m','LineWidth',1.5)
title('Overview')
xlabel('time locked to feedback [ms]')
ylabel('z-score normalized')
legend('ERP @ Cz','grand GFP','grand GMD')
grid on

%% 5. ERP image + smoothing
ch = 20; subject = 1;
gauss_std = 2;
data = squeeze(DATASET.data(ch,:,:,subject))';
data_nonerror = data(DATASET.labels(:,subject)==-1,:);
data_error = data(DATASET.labels(:,subject)==1,:);
% smooth with 2d gaussian kernel
data_ = imgaussfilt(data,gauss_std);
data_nonerror_ = imgaussfilt(data_nonerror,gauss_std);
data_error_ = imgaussfilt(data_error,gauss_std);

figure(3)
subplot(2,2,[1 3])
contourf(DATASET.times,1:150,data_,10,'linecolor','none')
xlabel('Time [ms]'); ylabel('trial')
title('ERP image channel Cz (all)')
colorbar

subplot(2,2,2)
contourf(DATASET.times,1:size(data_nonerror_,1),data_nonerror_,10,'linecolor','none')
xlabel('Time [ms]'); ylabel('trial')
title('ERP image channel Cz (non-error)')
colorbar
subplot(2,2,4)
contourf(DATASET.times,1:size(data_error_,1),data_error_,10,'linecolor','none')
xlabel('Time [ms]'); ylabel('trial')
title('ERP image channel Cz (error)')
colorbar


%% 6. Topographic Plots
timeslots = [226 277 351 609 871];
cmin = min(min(ERP.diff_Gavg)); cmax = max(max(ERP.diff_Gavg));
figure(4)
for t = 1:length(timeslots)
    subplot(1,length(timeslots),t)
    topoplot(ERP.diff_Gavg(:,find(floor(DATASET.times)==timeslots(t))),DATASET.chanlocs,'maplimits',[cmin cmax])
    title(strcat(num2str(timeslots(t)),' ms'))
end


%% 8. ROI-based statistical analysis
% N200: 180-325ms over frontocentral channels
% P300: 300-600ms over central channels

sel_chans = [2 6 10 11 19 20 24];
sel_times_n200 = DATASET.times>=180 & DATASET.times<=325;
sel_times_p300 = DATASET.times>=300 & DATASET.times<=600;

pool_n200 = squeeze(mean(mean(ERP.diff_avg(sel_chans,sel_times_n200,:),2),3));
pool_p300 = squeeze(mean(mean(ERP.diff_avg(sel_chans,sel_times_p300,:),2),3));

figure(5)
boxplot([pool_n200 pool_p300],{'N200','P300'})
ylabel('Amplitude [µV]')
title('ROI-based distributions (n=17)')
grid on

% statistical testing: t-test
[parametric.h_n200 parametric.p_n200] = ttest(pool_n200);
[parametric.h_p300 parametric.p_p300] = ttest(pool_p300);


%% 9. ERP-map and masking with coefficient of determination R²

% compute coefficient of determination for each subject
for ch=1:size(DATASET.data,1)
    for t=1:size(DATASET.data,2)
        for s=1:size(DATASET.data,4)
            r(ch,t,s) = corr(squeeze(DATASET.data(ch,t,:,s)),DATASET.labels(:,s));
        end
    end
    ch
end

R2 = r.^2;
R2_m = mean(R2,3);

subplot(1,3,1)
contourf(ERP.diff_Gavg(:,sel_time),20,'linecolor','none')
xticklabels = -200:200:1200;
xticks = linspace(1, size(ERP.diff_Gavg(:,sel_time), 2), numel(xticklabels));
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels)
set(gca, 'YTick', 1:27, 'YTickLabel', DATASET.chanlist)
caxis([-3 4])
ylabel('channel #')
xlabel('time locked to feedback [ms]')
title('spatiotemporal ERP map')
colorbar

subplot(1,3,2)
contourf(R2_m(:,sel_time),20,'linecolor','none')
xticklabels = -200:200:1200;
xticks = linspace(1, size(ERP.diff_Gavg(:,sel_time), 2), numel(xticklabels));
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels)
set(gca, 'YTick', 1:27, 'YTickLabel', DATASET.chanlist)
ylabel('channel #')
xlabel('time locked to feedback [ms]')
title('coeff. of determ. R²')
colorbar

subplot(1,3,3)
contourf(ERP.diff_Gavg(:,sel_time).*(R2_m(:,sel_time)>0.05),20,'linecolor','none')
xticklabels = -200:200:1200;
xticks = linspace(1, size(ERP.diff_Gavg(:,sel_time), 2), numel(xticklabels));
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels)
set(gca, 'YTick', 1:27, 'YTickLabel', DATASET.chanlist)
caxis([-3 4])
ylabel('channel #')
xlabel('time locked to feedback [ms]')
title('masked ERP map (R²>0.05)')
colorbar