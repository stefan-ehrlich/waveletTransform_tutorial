%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Tutorial II part C - sample solution
% Author: Stefan Ehrlich
% Date: 21.01.2020
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

addpath(genpath(strcat(pwd,'\functions\')));

inpath = strcat(pwd,'\data\');
load(strcat(inpath,'TF_DATA.mat'));

%% 1. Visualization of the TF data
% compute condition-wise subject averages
for s = 1:17
    nonerror_mean(:,:,:,s) = mean(TF_DATA.magnitude(:,:,:,TF_DATA.labels(:,s)==-1,s),4);
    error_mean(:,:,:,s) = mean(TF_DATA.magnitude(:,:,:,TF_DATA.labels(:,s)==1,s),4);
end
% compute grand average difference
diff_gAvg = mean(error_mean,4)-mean(nonerror_mean,4);

% visualization
tfviewerx(TF_DATA.times,TF_DATA.frex,double(diff_gAvg),TF_DATA.chanlocs)

%% 2. Univariate statistical testing without correction for multiple comparisons

for ch=1:size(nonerror_mean,1)
    for f = 1:size(nonerror_mean,2)
        for t = 1:size(nonerror_mean,3)
            [h(ch,f,t) p(ch,f,t)] = ttest(squeeze(nonerror_mean(ch,f,t,:)),squeeze(error_mean(ch,f,t,:)));
        end
    end
    ch
end

statsmask_005 = p<0.05; % mask for p-level threshold of 0.05
statsmask_0001 = p<0.001; % mask for p-level threshold of 0.001

tfviewerx(TF_DATA.times,TF_DATA.frex,double(diff_gAvg).*statsmask_0001,TF_DATA.chanlocs)

%% 3. Univariate statistical testing with Bonferroni correction
nComp = size(nonerror_mean,1)*size(nonerror_mean,2)*size(nonerror_mean,3);
p_bonferroni = 0.05/nComp;

statsmask_bonferroni = p<p_bonferroni;

tfviewerx(TF_DATA.times,TF_DATA.frex,double(diff_gAvg).*statsmask_bonferroni,TF_DATA.chanlocs)

%% Non-parametric permutation testing with correction for multiple comparisons
nPerms = 250;
tic;
for i=1:nPerms
    % step 1: shuffle the labels
    for s=1:17
        shuffleMat = randperm(150);
        labels_shuffled(:,s) = TF_DATA.labels(shuffleMat,s);
    end

    % step 2: compute surrogate condition wise averages
    for s = 1:17
        nonerror_mean_shuffled(:,:,:,s) = mean(TF_DATA.magnitude(:,:,:,labels_shuffled(:,s)==-1,s),4);
        error_mean_shuffled(:,:,:,s) = mean(TF_DATA.magnitude(:,:,:,labels_shuffled(:,s)==1,s),4);
    end
    % step 3: compute surrogate grand average difference
    diff_gAvg_shuffled = mean(error_mean_shuffled,4)-mean(nonerror_mean_shuffled,4);

    % step 4: store max and min values
    maxVal(i) = max(max(max(diff_gAvg_shuffled)));
    minVal(i) = min(min(min(diff_gAvg_shuffled)));
    i
end
elapsed=toc;

% step 5: compute global surragate min and max values
pctmin = prctile(minVal,2.5)
pctmax = prctile(maxVal,97.5)
figure(1)
subplot(2,1,1)
hist(maxVal)
hold on
plot([pctmax pctmax],[0,nPerms/5],'g','LineWidth',1.5)
title('max-value distribution')
subplot(2,1,2)
hist(minVal)    
hold on
plot([pctmin pctmin],[0,nPerms/5],'g','LineWidth',1.5)
title('min-value distribution')

% step 6: comparing original grand average difference with global surrage min and max
% generate stats mask
statsmask_perm = (diff_gAvg>pctmax | diff_gAvg<pctmin);

tfviewerx(TF_DATA.times,TF_DATA.frex,double(diff_gAvg).*statsmask_perm,TF_DATA.chanlocs)



