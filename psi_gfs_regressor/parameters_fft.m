function [ X, YO, YO_Im ] = parameters_fft(x, Fs, freq_band, plt)

Np = length(x);

% FFT analysis

NFFT = 2 ^ nextpow2(Np); % Next power of 2 from length of data
YO_Im = fft(x, NFFT) ./ Np; YO_Im = YO_Im(1:NFFT/2 + 1); 
YO = 2 * abs(YO_Im); YO = YO';

% Axis

fmin = freq_band(1); fmax = freq_band(2);
X = Fs / 2 * linspace(0, 1, NFFT/2 + 1); 
X_aux = logical( ( X <= fmax ) .* ( X >= fmin ) ); X = ones(1,1) * X(X_aux);

YO = YO(logical(X_aux)); YO_Im = YO_Im(logical(X_aux));

if plt
    % No averaging
    figure('Name', 'PSD'), plot(X, YO, 'color', 0.5 * rand * [ 1 1 1 ], 'linewidth', 2); hold on; grid on;
    set(gca, 'YLim', [ 0 4 ]); xlabel('Frequency [Hz]'); ylabel('Power');
end