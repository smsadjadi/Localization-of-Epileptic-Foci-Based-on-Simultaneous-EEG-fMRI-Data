close all; clear; clc;

%% IED template extraction

SubDataPath = 'D:\Academic\Datasets\EpilepsyEEGfMRI\Sub_01\';
outsideEEG = pop_loadset([SubDataPath,'EEG\Processed\Final\outside_preprocessed_IEDs.set']);
data = outsideEEG.data([1:31,33:64],:,:);
nbchan = outsideEEG.nbchan-1; % One channel is belongs to ECG.
trials = outsideEEG.trials; % Number of events
IED = zeros(nbchan,250,trials); % The IED takes 1 sec.
t = linspace(-.5,.5,250);
spike = zeros(nbchan,250);

% Channel Spike Template
spikeChannel = [];
for tr = 1:trials
    signal = data(:,:,tr);
    sigcenter = signal(:,231:270);
    [chpeak,tpeak] = find(abs(sigcenter) == max(abs(sigcenter(:))));
    tpeak = tpeak+230;
    IED(:,:,tr) = signal(:,tpeak-125:tpeak+124);
    spike = spike+IED(:,:,tr);
    spikeChannel = [spikeChannel,outsideEEG.chanlocs(chpeak).labels,','];
    subplot(2,3,[2,3]); plot(t,IED(:,:,tr)',':'); hold on;
end

spike = spike/trials;
spikePattern = mean(spike);
subplot(2,3,[2,3]); plot(t,spikePattern,'y','LineWidth',2);
subplot(2,3,[5,6]); plot(t,spikePattern);
subplot(2,3,[1,4]); eegplot(mean(IED,3),250,0,0.0002,0.321)

% % Channel Spike template
% spikeChannel = []; i = 0;
% for tr = 1:trials
%     signal = data(:,:,tr);
%     [chpeak,tpeak] = find(signal == max(signal(:)));
%      if tpeak>150 && tpeak<350
%         i = i+1;
%         IED(:,:,tr) = signal(:,tpeak-125:tpeak+125-1);
%         spike = spike+IED(:,:,tr);
%         spikeChannel = [spikeChannel,outsideEEG.chanlocs(chpeak).labels,','];
%      end
% end
% spike = spike/i;
% plot(mean(spike));

%% IED detection from inside the MR scanner

insideEEG = pop_loadset([SubDataPath,'EEG\Processed\Final\inside_fMRIb_BCG_ICA.set']);
data = insideEEG.data;

% Cross Correlation
cc = zeros(nbchan,size(spike,2)+size(data,2));
for ch = 1:nbchan
    c = xcorr(spike(ch,:),data(ch,:));
    cc(ch,:) = c(1:size(spike,2)+size(data,2));
end

% Plot
ccSorted = sort(cc,'descend');
ccSum = sum(ccSorted(1:10,:));
figure; plot(linspace(0,size(cc,2)/250,size(cc,2)),ccSum);

% Thresholding
% hist(ccSum);
ccSumSorted = sort(ccSum);
thresh = ccSumSorted(end-floor(0.005*length(ccSumSorted)));
ccSumTh = ccSum;
ccSumTh(ccSumTh<thresh) = 0; ccSumTh(ccSumTh~=0) = 1;
figure; plot(linspace(0,size(cc,2)/250,size(cc,2)),ccSumTh); ylim([-0.5,1.5]);

% Generating Spike series
SpikeTiming = ccSumTh;
save([SubDataPath,'Regressors\SpikeTiming.mat'],'SpikeTiming');

%% Generating linear regressor

TR = 2.5;
samplingRate = 250;
fMRIvolumes = 117;
DeleteVolumes = 9;
IEDtimes = find(SpikeTiming==1)/samplingRate;
IEDtimesDiff = IEDtimes - [0,IEDtimes(1:end-1)];
IEDtimes(IEDtimesDiff<1) = [];
% IEDtimes(IEDtimes>fMRIvolumes*TR) = [];
% IEDtimes(IEDtimes<DeleteVolumes*TR) = [];
IEDduration = 0.1*ones(size(IEDtimes));
IEDstate = 1*ones(size(IEDtimes));
LinearRegressor = [IEDtimes;IEDduration;IEDstate]';
dlmwrite([SubDataPath,'Regressors\LinearRegressor.txt'],LinearRegressor,'delimiter',' ')

%% Alpha base extraction

outsideEEG = pop_loadset('D:\Academic\Datasets\EpilepsyEEGfMRI\Sub_01\EEG\Processed\outside\outside-rs250-LPfilter-cleanline&notch-ICA-reduce.set');
data = outsideEEG.data([1:31,33:64],:,:);
times = outsideEEG.times;
nbchan = outsideEEG.nbchan-1;
trials = outsideEEG.trials;
p = data.^2;
p = sort(p,'descend');
pAlpha = sum(p(1:10,:));
plot(linspace(0,size(pAlpha,2)/250,size(pAlpha,2)),pAlpha);
% IED Base Detection
% In the alpha base, there is a disturbance with the period of 0.1 sec in each side of IED spike.
pAlphaTh = zeros(size(pAlpha));
for i=251:length(pAlpha)-250
    if sum(pAlpha(i-25:i+24))>sum(pAlpha(i-75:i-26)) && ...
       sum(pAlpha(i-25:i+24))>2*sum(pAlpha(i-150:i-76)) && ...
       sum(pAlpha(i-25:i+24))>sum(pAlpha(i+25:i+74)) && ...
       sum(pAlpha(i-25:i+24))>2*sum(pAlpha(i+75:i+149))
       pAlphaTh(i) = 1;
    end
end
% plot(linspace(0,size(pAlpha,2)/250,size(pAlpha,2)),pAlphaTh); ylim([-1,3]);
len = min(length(ccSumTh),length(pAlphaTh));
figure; plot(linspace(0,len/250,len),ccSumTh(1:len)); ylim([-1,3]);
hold on; plot(linspace(0,len/250,len),pAlphaTh(1:len)); ylim([-1,3]); legend('Spike','Base');
%Thresholding Alpha wave power
pAlphaSorted = sort(pAlpha);
pthresh = pAlphaSorted(end-floor(0.01*length(pAlphaSorted)));
pAlphaTh = pAlpha;
pAlphaTh(pAlphaTh<pthresh) = 0; pAlphaTh(pAlphaTh~=0) = 1;
plot(linspace(0,size(pAlpha,2)/250,size(pAlpha,2)),pAlphaTh); ylim([-1,3]);
len = min(length(ccSumTh),length(pAlphaTh));
figure; plot(linspace(0,len/250,len),ccSumTh(1:len)); ylim([-1,3]);
hold on; plot(linspace(0,len/250,len),pAlphaTh(1:len)); ylim([-1,3]);