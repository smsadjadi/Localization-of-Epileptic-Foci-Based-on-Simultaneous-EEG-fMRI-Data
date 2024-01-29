function eegplot(signal,fs,t0,Vscale,yscale)
    % signal = mean(IED,3); % signal = IED(:,:,2); fs = 250; t0 = -0.5; Vscale = 0.0002; yscale = 0.321;
    nbchan = size(signal,1); nbsamp = size(signal,2);
    t = linspace(t0, t0+nbsamp/fs, nbsamp);                      % Time Vector
    EEG = signal*Vscale;                                         % Simulate EEG (mV)
    ofst = (1:nbchan)*0.005;                                     % ‘Offset’ Vector
    EEGp = bsxfun(@plus, EEG', ofst)';                           % Add ‘Offset’ To Each Row
    plot(t, EEGp,'color','b')                                    % Plot EEG
    axis([xlim 0 yscale])                                        % Set Axis Limits
    ChC = regexp(sprintf('Ch-%02d ',1:nbchan),' ','split');      % Y-Tick Labels
    yt = ofst+0.0025;                                            % Y-Yick Positions
    set(gca,'YTick',yt,'YTickLabel',ChC(end-1:-1:1))             % Set Tick Labels
end