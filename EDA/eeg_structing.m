clc; clear;

eeglabData = pop_loadset('D:\Academic\Datasets\EpilepsyEEGfMRI\Sub_01\EEG\Processed\inside\Inside-5min(1)-fmrib-rs250-cleanline-QRSremove-HP1Hz_chLoc.set');
EEGchan = [1:31,33:64];
EEG = struct();
EEG.data = double(eeglabData.data(EEGchan,:));
EEG.labels = {eeglabData.chanlocs(EEGchan).labels}';
EEG.type = 'EEG';
EEG.nbchan = eeglabData.nbchan;
EEG.points = length(eeglabData.times);
EEG.srate = eeglabData.srate;
EEG.labeltype = 'standard';
EEG.unit = 'uV';

sampleLength = [1:2500];
EEG.data = EEG.data(:,sampleLength);
EEG.points = length(sampleLength);
EEG.nbchan = size(EEG.data,1);

save('EEGdata.mat','EEG','-mat');



