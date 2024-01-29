function [reg_fit, reg_fit_td, reg_fit_dd] = fit2reg(data, trigs, fs, NTR)
    
    % fit2reg() - Conversion of fitting results into fMRI regressor
    % 
    % INPUTS:
    %   data - results from the fitting (should be a five column matrix)
    %   trigs - MR volumes triggers in time frame
    %   fs - EEG sampling rate in Hertz
    %   NTR - number of volume
    %
    % OUTPUTS:
    %   reg_fit - regressor obtained after convolution by the HRF
    %   reg_fit_td - time derivative
    %   reg_fit_dd - dispersion derivative
    
    trigs(end+1)=trigs(end)+round(mean(diff(trigs)));
    
    fit_corr = data(:,5);
    mean_r=mean(fit_corr(trigs(1):trigs(end)))
    std_t=std(fit_corr(trigs(1):trigs(end)))
    fit_corr = fit_corr.^2;
    hrf= spm_hrf(1/fs);
    hrf=hrf/std(hrf);
    td = diff(hrf);
    td=td/std(td);
    dd = diff(td);
    dd=dd/std(dd);
    
    reg = conv(fit_corr,hrf);
    reg_td = conv(fit_corr,td);
    reg_dd = conv(fit_corr,dd);
    
    reg_fit = zeros(NTR,1);
    reg_fit_td = zeros(NTR,1);
    reg_fit_dd = zeros(NTR,1);
    
    
    for i=1:NTR
        reg_fit(i)=mean(reg(trigs(i):trigs(i+1)-1));
        reg_fit_td(i)=mean(reg_td(trigs(i):trigs(i+1)-1));
        reg_fit_dd(i)=mean(reg_dd(trigs(i):trigs(i+1)-1));
    end
    
    figure
    plot(reg_fit)

end