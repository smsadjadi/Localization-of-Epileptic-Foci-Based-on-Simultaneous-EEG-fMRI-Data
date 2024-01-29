function [ vf, wavelet_cf, power ] = tf_decomposition_linear_scale(fsample, data, ic_or_channel, filters, plot_TF)
% EEG -> data(n_channels * n_frames)
% ic_or_channel -> index of the Independent Component OR Channel in which
% the Time-Frequency decomposition will be applied
% is_ic -> binary: 1-ic, 0-channel

R = 7;% 3600; % Morlet Wavelet factor
maxf = filters(2); % 93.74; %impF(end); % S.fsample / 2 or other;
minf = filters(1); % 93.32; % Defined by me
B = 50; % 25; # Frequency Samples
% fsample = dataset.srate;

AD = data(ic_or_channel, :);
vf = linspace(minf, maxf, B); % array containing the frequency bins (more frequency bins for the lowest frquencies)

dataeegF = fft(AD);
N = length(dataeegF);
wavelet_cf = zeros(length(vf), N); % morlet wavelet coefficients
n = 1;

for F = vf
    Fwavelet = single(MorletWavelet(F, R, fsample, N)); % compute the Morlet wavelet at a given frequency bin F
    wavelet_cf(n, :) = ifft(dataeegF .* fft(Fwavelet)) ./ N; %convolution in the frequency domain
    n = n + 1;
end

power = abs(wavelet_cf) .^ 2; % power at a given time t around a frequency f -> P(f,t)

if plot_TF
    % Plot Power Spectrum w/ TF + data
    Dtf_f = abs(wavelet_cf) + .001; % Mexer aqui para alterar escala de cores
    
    figure('Name', 'TF', 'Position', [10, 550, 1900, 400], 'NumberTitle','off')
    subplot(211)
    imagesc((1:size(power,2)) ./ fsample, 1:length(vf), squeeze(Dtf_f));
    Fplot = [ 1, 5:5:45 ];
    hold on
    plot([1 N],[Fplot', Fplot'],'k')
    hold off
    set(gca,'YTick',Fplot)
    set(gca,'YTickLabel', [ 1, 5:5:45 ], 'FontSize', 10)
    ylabel('Frequency (Hz)','FontSize', 10);
    cbar;
    subplot(212)
    plot((1:size(power,2)) ./ fsample, AD)
    axis('tight')
    xlabel('Time [s]', 'FontSize', 10) 
end

return;