close all; clear; clc;

%% IED template extraction

SubDataPath = 'D:\Academic\Datasets\EpilepsyEEGfMRI\Sub_01\';
outsideEEG = pop_loadset([SubDataPath,'EEG\Processed\Final\outside_preprocessed_IEDs.set']);
data = outsideEEG.data([1:31,33:64],:,:);
nbchan = outsideEEG.nbchan-1; % One channel is belongs to ECG.
trials = outsideEEG.trials-1; % Number of events %%%%%%%%%%%%%%%%%%%%% -1 %%%%%%%%%%%%%%%%%%%%%
IED = zeros(nbchan,250,trials); % The IED takes 1 sec.
spike = zeros(nbchan,250);

% Spike pattern
spikeChannel = [];
for tr = 1:trials
    signal = data(:,:,tr);
    sigcenter = signal(:,231:270);
    [chpeak,tpeak] = find(abs(sigcenter) == max(abs(sigcenter(:))));
    tpeak = tpeak+230;
    IED(:,:,tr) = signal(:,tpeak-125:tpeak+124);
    spike = spike+IED(:,:,tr);
    spikeChannel = [spikeChannel,outsideEEG.chanlocs(chpeak).labels,','];
    subplot(2,3,[2,3]); plot(IED(:,:,tr)',':'); hold on;
end

spike = spike/trials;
spikePattern = mean(spike);
subplot(2,3,[2,3]); plot(spikePattern,'y','LineWidth',2);
subplot(2,3,[5,6]); plot(spikePattern);
subplot(2,3,[1,4]); eegplot(mean(IED,3),250,0,0.0002,0.321)

% % Spike pattern
% spikeChannel = []; i = 0;
% times = outsideEEG.times; % Time period (msec)
% for tr = 1:trials
%     for ch = 1:nbchan
%         signal = data(ch,:,tr);
%         [chpeak,tpeak] = find(signal == max(signal));
%         if tpeak>200 && tpeak<300
%             i = i+1;
%             IED(ch,:,tr) = signal(tpeak-125:tpeak+125-1);
%             subplot(211); plot(times(size(times,2)/2-125:size(times,2)/2+125-1),IED(ch,:,tr)); hold on
%             spikePattern = spikePattern+IED(ch,:,tr);
%             spikeChannel = [spikeChannel,outsideEEG.chanlocs(chpeak).labels,','];
%         end
%     end
% end
% spikePattern = spikePattern/i;
% subplot(211); plot(times(size(times,2)/2-125:size(times,2)/2+125-1),spikePattern,'*'); xlabel('time (ms)'); ylabel('voltage (uV)');
% subplot(212); plot(times(size(times,2)/2-125:size(times,2)/2+125-1),spikePattern); xlabel('time (ms)'); ylabel('voltage (uV)');

% Spike timing
% load([SubDataPath,'Regressors\SpikeTiming.mat']);

%% eIC detection

insideEEG = pop_loadset([SubDataPath,'EEG\Processed\Final\inside_fMRIb_BCG_ICA.set']);
data = insideEEG.data;
components = insideEEG.icaact;

if isempty(components)
    insideEEG.icaact = (insideEEG.icaweights*insideEEG.icasphere)*insideEEG.data(insideEEG.icachansind,:);
    components = insideEEG.icaact;
end

% eIC detection (method 1)
eIC_num = 1:10;
disp(['eIC number = ',num2str(eIC_num)])

% % eIC detection (method 2)
% cc = zeros(size(components,1),size(spike,2)+size(components,2));
% CC = zeros(1,size(components,1));
% for ch = 1:size(components,1)
%     c = xcorr(spike(ch,:),components(ch,:));
%     cc(ch,:) = c(1:size(spike,2)+size(components,2));
%     CC(ch) = sum(cc(ch,:).^2);
% end
% eIC_num = find(CC == max(CC));
% disp(['eIC number = ',num2str(eIC_num)])
% 
% % eIC detection (method 3)
% load([SubDataPath,'Regressors\SpikeTiming.mat']);
% SpikeSeries = conv(SpikeTiming,spikePattern);
% SpikeSeries = SpikeSeries(floor(linspace(1,length(SpikeSeries),length(components))));
% figure; plot(SpikeSeries)
% for ch = 1:size(components,1)
%     c = xcorr(SpikeSeries,components(ch,:));
%     CC(ch) = sum(c.^2);
% end
% eIC_num = find(CC == max(CC));
% disp(['eIC number = ',num2str(eIC_num)])

%% Generating component regressor

TR = 2.5;
samplingRate = 250;
fMRIvolumes = 117;
DeleteVolumes = 9;
for i=1:length(eIC_num)
    eIC = components(eIC_num(i),:);
   figure; plot(eIC);
    ICregressor = zeros(1,fMRIvolumes-DeleteVolumes);
    tr = floor(linspace(1,length(eIC),fMRIvolumes-DeleteVolumes+1));
    for ii = 1:length(tr)-1
        onetr = eIC(tr(ii):tr(ii+1)-1);
        [e,t] = max(abs(onetr));
        ICregressor(ii) = onetr(t);
    end
    dlmwrite([SubDataPath,'Regressors\Component',num2str(eIC_num(i)),'Regressor.txt'],ICregressor')
   figure; plot(ICregressor);
end

% for i=1:length(eIC_num)
%     eIC = components(eIC_num(i),:);
%     ICregressor = eIC(floor(linspace(1,length(eIC),fMRIvolumes-DeleteVolumes)));
%     dlmwrite([SubDataPath,'Regressors\Component',num2str(eIC_num(i)),'Regressor.txt'],ICregressor')
% end