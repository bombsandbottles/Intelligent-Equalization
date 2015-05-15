% DSTII Final Project %
% Harrison Zafrin hzz200
% filtered_spectrum = output smoothed magnitudes spectrum
% size = size of the moving average filter
% X_mag_avg = avg spectrum of all mag spectrums together
% -------------------------------------------------------------------------
% Manual homebrew implementation of centerd moving average filter
% -------------------------------------------------------------------------
function [ filtered_spectrum ] = movingavgfilter_17pnt( X_mag_avg, size, fftparams, fs)

temp = X_mag_avg;

% -------------------------------------------------------------------------
% This section takes the moving average of the entire signal as normal
% -------------------------------------------------------------------------

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

% -------------------------------------------------------------------------
% This section then appends the original <200hz values back and keeps the 
% >200hz values smoothed
% -------------------------------------------------------------------------

% Get win_size to get hertz resolution of FFT
win_size = fftparams.win_size;

% Create hertz vector
hertz_vals = linspace(0, fs, win_size);

% Find the index value where X_mag_avg(k) is above 200hz
k = 1;
while hertz_vals(k) < 200
    k = k+1;
end

% Values of X_mag_avg below 200hz
X_b200 = temp(1:k-1);

% Values of X_mag_avg above 200hz
X_a200 = filtered_spectrum(k:end);

filtered_spectrum = [X_b200' X_a200];

end