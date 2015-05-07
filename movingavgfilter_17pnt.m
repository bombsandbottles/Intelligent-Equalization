% DSTII Final Project %
% Harrison Zafrin hzz200
% X_mag_avg = avg spectrum of all mag spectrums together
% -------------------------------------------------------------------------
% Manual implementation of centerd moving average filter
% -------------------------------------------------------------------------
function [ filtered_spectrum ] = movingavgfilter_17pnt( X_mag_avg, size)

% Pad amount
pad_amt = (size-1)/2;

% Get zeros to append for y-1->y-8
left = zeros(1, pad_amt);
right = zeros(1, pad_amt);

% Concat Them All to zero-pad front and back
X_mag_avg = [left X_mag_avg' right];

% Offline buffer for moving average based on size
win_size = size;
hop_size = 1;

% Get amount of sample overlap per window
n_overlap = win_size - hop_size;

% Buffer with n_overlap
X_mag_avg = buffer(X_mag_avg, win_size, n_overlap, 'nodelay');

% Use the Power of MATLAB to avg the windows
filtered_spectrum = mean(X_mag_avg);

end