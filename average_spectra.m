% DSTII Final Project %
% Harrison Zafrin hzz200
% x_t = time domain signal
% X_mag_mean = average magnitude response of x_t
% X_mag_norm = normalized magnitude response so that the bin sum is 1
% X_mag_cum = cumulative distribution function of X_norm
% fs = sampling rate
% -------------------------------------------------------------------------
% This function takes a 4096 point STFT of x_t and returns the averaged magnitude
% spectra across time into 4096 bins with about 10hz resolution per bin
% -------------------------------------------------------------------------
function [ X_mag, X_mag_mean, X_mag_cum, fs ] = average_spectra( filepath, fftparams, discrete_signal )

% Perform averaging on a .wav
if nargin < 3
    % Import Audio
    [x_t, unused, fs, t] = import_audio_stereo( filepath );

    % lol
    unused = unused/2;
    
% Perform averaging on a currently discretized signal for testing
elseif nargin == 3
    
    filepath = 0;
    fs = 44100;
    x_t = discrete_signal;
    
end

% Find the length of x_t
x_len = length(x_t);

% -------------------------------------------------------------------------
% STFT
% -------------------------------------------------------------------------

% Define STFT Params
win_size = fftparams.win_size;
hop_size = win_size/2;

% Get amount of sample overlap per window
n_overlap = win_size - hop_size;

% Create Window a hanning window to prevent spectral leakage
window = hann(win_size);

% Buffer x_t with n_overlap
x_t_buff = buffer(x_t, win_size, n_overlap, 'nodelay');

% Create Window Matrix
window_mat = repmat(window, 1, size(x_t_buff, 2));

% Window the Signal
x_t_windowed = x_t_buff .* window_mat;

% Get the magnitude response via STFT and subsequent removing of phase
X_mag = abs(fft(x_t_windowed));

% Remove Mirror Image past fs/2
X_mag = X_mag(1:end/2, :);

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
% Normalization of Spectra
% -------------------------------------------------------------------------

% Scale so bin sum would be 1
X_mag_norm = X_mag_mean/sum(X_mag_mean);

% Accumulate Over The Bins (cumulative distribution function)
X_mag_cum = cumsum(X_mag_norm);
X_mag_cum = X_mag_cum';

end

% -------------------------------------------------------------------------
% Test Code IGNORE
% -------------------------------------------------------------------------
% 
% % Create the Average Spectra For Multiple Normalized Spectra
% song_num = 1;
% 
% % Mean calculation of each point in the cumulative distribution
% X_cum_mean = mean(X_mag_cum', 2);
% 
% % Create a Vector of Xc(k) - Xc(k-1)
% for k=2:length(X_cum_mean)
%     
%     % Compute the difference between adjacent bins    
%     X_cum_mean(k) = X_cum_mean(k) - X_cum_mean(k-1);
%     
% end
% 
% % Compute the average spectra across multiple spectra
% X_avg = (sum(X_mag_mean, 2)/1) .* X_cum_mean;

