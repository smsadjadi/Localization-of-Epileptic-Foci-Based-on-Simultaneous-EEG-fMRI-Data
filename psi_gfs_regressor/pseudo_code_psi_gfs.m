% pseudo-code

% --------------- PSI --------------- %

% FIRST STEP: 
% - select the channel with the highest amplitude of the average IED
% - store the index of that channel as ch

% SECOND STEP: compute the time-frequency decomposition of that reference channel
% fsample: sampling frequency

% data_sync: EEG data: N x M matrix (N - number of channels, M - number of time points)
%            data_sync has been already re-referenced using a Laplacian (very important step)

% ch: index of the reference channel

% filters: vector [ f1, f2 ], indicating the range of frequencies [in Hz] in the EEG data after the band-pass filtering
%          if filters = [ 1, 45 ], then the EEG data was band-pass filtered between 1 and 45 Hz

% vf: frequency bins on which the TF decomposition was computed
% wt_ref: Morlet wavelet coefficients (complex numbers)

[ vf, wt_ref, ~ ] = tf_decomposition(fsample, data_sync, ch, filters, 0);

% get the indices of vf spanning the frequency of interest for PSI
% freq_band: vector [ f1, f2 ]; if freq_band = [ 1 30 ], the average PSI will be computed across the frequency bins spanning the 1-30 Hz range
fi = find(vf <= freq_band(1), 1, 'last');
fe = find(vf <= freq_band(2), 1, 'last'); if isempty(fi), fi = 1; end

% initialize some variables:
% nb_ch = number of channels
% S = number of fMRI time points
% F = length(vf) [number of frequency bins]
% N = number of EEG samples in a single TR (N = fsample * TR, with TR in seconds)
ph_ref = angle(wt_ref); F = length(vf); W = zeros(nb_ch, S, F);
PSI = zeros(nb_ch - 1, S, F); m_PSI = zeros(nb_ch - 1, S); max_stdPSI = zeros(nb_ch - 1, S);

% THIRD STEP: compute the PSI between the reference channel and all the other channels 
% for segments equal to the length of each TR -> will result in one PSI value per TR, ready to be used in a GLM of the fMRI data
    
count = 0;
for n2 = 1:nb_ch
    if n2 ~= ch
        
        count = count + 1;
        [ ~, wt, ~ ] = tf_decomposition_linear_scale(dataset_bcg, data_sync, n2, filters, 0);
        ph = angle(wt); 
        
        % get the cross wavelet coefficients between the two channels
        xwt = wt_ref .* conj(wt); wg = abs(imag(xwt));
        
        for j = 1:S
            
            % compute the phase difference on the segment j
            DPH = ph_ref(:, (j - 1) * N + 1:j * N) - ph(:, (j - 1) * N + 1:j * N);
            
            % get the Y and X variables
            Y = (1 / N) .* sum(sin(DPH), 2); X = (1 / N) .* sum(cos(DPH), 2);
            
            % compute the PSI
            PSI(count, j, :) = sqrt(X .^ 2 + Y .^ 2);
            
            % get the average PSI across a frequency of interest
            m_PSI(count, j) = mean(squeeze(PSI(count, j, fi:fe)));
        end
        
        % get the combination of channels n2 and ch yielding the PSI with the highest variance
        [ ~, m ] = max(std(squeeze(PSI(count, :, fi:fe)), 0, 1));
        
        % store that PSI values
        max_stdPSI(count, :) = squeeze(PSI(count, :, m));
        
        % some plots
        if plt
            figure('Name', 'PSI')
            imagesc((1:S) .* dataset_bcg.TR_sec, 1:F, squeeze(PSI(count, :, :))');
            alfa = (filters(2) / filters(1)) ^ (1 / F) - 1;
            Fplot = (log([1 2 5 8 13 20 30 45 50 70 90]) - log(vf(1))) ./ log(1 + alfa) + 2;
            hold on, plot([ 1 N ], [ Fplot', Fplot' ], 'k'), hold off
            set(gca,'YTick',Fplot)
            set(gca,'YTickLabel', [1 2 5 8 13 20 30 45 50 70 90], 'FontSize', 10)
            ylabel('Frequency (Hz)','FontSize', 10);
        end
    end
end

% --------------- GFS --------------- %

% FFT decomposition of a single EEG segment of length N
[ X, ~, Yaux ] = parameters_fft(data_sync(1, 1:N), fsample, filters, 0);

% initialize some variables (F is the number of frequency bins)
F = length(Yaux); SS = zeros(nb_ch, S, F); CC = zeros(nb_ch, S, F); YY = zeros(nb_ch, S, F); GFS = zeros(S, F);

for j = 1:S
    
    % compute the FFT for each segment j (of length N) and for each channel n
    
    for n = 1:nb_ch
        
        % demean the data
        y = data_sync(n, :) - mean(data_sync(n, :));
        
        % get the FFT
        [ ~, ~, Y ] = parameters_fft(y((j - 1) * N + 1:j * N), dataset_bcg.srate, filters, 0);
        ph = angle(Y); SS(n, j, :) = sin(ph); CC(n, j, :) = cos(ph); YY(n, j, :) = Y;
        
    end
    
    % get the GFS for each frequency bin f, and each segment j
    
    for f = 1:F
        
        % demean (for pca)
        yy = squeeze(YY(:, j, f)) - mean(squeeze(YY(:, j, f)));
        
        % combine the real and imaginary parts of the FFT
        re = real(yy); im = imag(yy); vec = [ re'; im' ];
        
        % apply PCA with only two components, and get their eigen values
        [ ~, mapping ] = pca(vec', 2); eigval = mapping.lambda;
        
        % compute GFS
        GFS(j, f) = abs(eigval(1) - eigval(2)) / (eigval(1) + eigval(2));
    end
end

% get the average GFS across the frequency of interest (same rationale as for PSI, see above)
m_gfs = mean(GFS(:, fi:fe), 2);

% some plots
if plt
    figure('Name', 'GFS')
    imagesc((1:S) .* dataset_bcg.TR_sec, X, GFS);
    xlabel('Time [s]', 'FontSize', 10);
    ylabel('Frequency [Hz]','FontSize', 10); axis square
end

% convolution of GFS (and PSI) with HRF to then analyze fMRI data

% Building the HRF (Double Gamma) function with SPM
% TR is in seconds
HRF_EEG.dt = TR;
HRF_EEG.name = 'hrf';
HRF_EEG = spm_get_bf(HRF_EEG);

U.u = m_gfs; U.name = {'metrics'};
eeg_reg = spm_Volterra(U, HRF_EEG.bf)'; % convolution between HRF and EEG metric