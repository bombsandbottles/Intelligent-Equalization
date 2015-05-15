% DSTII Final Project %
% Harrison Zafrin hzz200
% x = A windowed time domain signal
% x_t = time domain vector version of x
% -------------------------------------------------------------------------
% Averages a spectra when given a windowed time domain matrix
% -------------------------------------------------------------------------
function [ averaged_spectra ] = average_spectra_matrix( x, x_t )

% Get win_size
win_size = size(x, 1);

% Get the magnitude response via STFT and subsequent removing of phase
X_mag = abs(fft(x));

% Remove Mirror Image past fs/2
X_mag = X_mag(1:end/2, :);

% Find the length of x_t
x_len = length(x_t);

% -------------------------------------------------------------------------
% Create Mean Magnitude over Time
% -------------------------------------------------------------------------

% Create numerator such that its the Sum of Xmag(k, t) over t
numerator = sum(X_mag, 2);

% Create denominator such that its x_len/win_size
if mod(x_len, win_size) == 0
    denominator = x_len/win_size;
else
    denominator = (x_len/win_size) + 1;
end

X_mag_mean = numerator/denominator;

% -------------------------------------------------------------------------
% Convert 2 dB
% -------------------------------------------------------------------------

% Scale so bin sum would be 1
X_mag_norm = X_mag_mean/sum(X_mag_mean);
% Accumulate Over The Bins (cumulative distribution function)
X_mag_cum = cumsum(X_mag_norm);
X_mag_cum = X_mag_cum';

[ X_mag_avg, freq_vector ] = create_overall_average_spectrum( X_mag_mean,...
                                                              X_mag_cum,...
                                                              size(X_mag_cum, 1),...
                                                              44100);

% Scale to dB
averaged_spectra = mag2db(X_mag_avg);

end

