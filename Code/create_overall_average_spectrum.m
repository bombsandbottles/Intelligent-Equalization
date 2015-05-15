% DSTII Final Project %
% Harrison Zafrin hzz200
% mag_mean_matrix = matrix containing the individually averaged spectrums
% mag_cum_matrix = cumulative distribution function matrix of magnitudes
% song_num = number of songs to average
% X_mag_avg = avg spectrum of all mag spectrums together
% fs = sampling rate
% -------------------------------------------------------------------------
% Compute Average Spectrum from normalized cumulative distributions
% -------------------------------------------------------------------------
function [ X_mag_avg, freq_vector ] = create_overall_average_spectrum( mag_mean_matrix,...
                                                                       mag_cum_matrix,...
                                                                       song_num,...
                                                                       fs)
% Mean calculation of each point in the cumulative distribution
mag_cum_matrix = mean(mag_cum_matrix', 2);

% Allocate k-1 at first
prev_mag = 0;

% Create a Vector of Xc(k) - Xc(k-1)
for k=1:length(mag_cum_matrix)
    
    % Compute the difference between adjacent bins    
    temp(k) = mag_cum_matrix(k) - prev_mag;
    
    % Update so that previous mag is now mag_cum_matrix(k-1)
    prev_mag = mag_cum_matrix(k);
    
end

% Compute the average spectra across multiple spectra
X_mag_avg = (sum(mag_mean_matrix, 2)/song_num) .* temp';

% % Create Frequency Vector
freq_vector = linspace(0, fs/2, length(X_mag_avg));

end

