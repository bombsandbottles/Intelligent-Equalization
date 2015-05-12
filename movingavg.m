% DSTII Final Project %
% Harrison Zafrin hzz200
% X_mag_avg = avg spectrum of all mag spectrums together
% -------------------------------------------------------------------------
% Manual implementation of centerd moving average filter, only above 200hz
% -------------------------------------------------------------------------
function [ filtered_spectrum ] = movingavg( X_mag_avg, size, fftparams, fs)

temp = X_mag_avg;

% -------------------------------------------------------------------------
% This section takes the moving average of the entire signal as normal
% -------------------------------------------------------------------------

% Create mavg params
b = (1/size)*ones(1,size);
a = 1;

% Filter those vals above 200hz
y = filter(b, a, X_mag_avg);

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
X_a200 = y(k:end);

% Create mavg params
b = (1/size)*ones(1,size);
a = 1;

% Put them back together
filtered_spectrum = [X_b200 ; y];

end


% -------------------------------------------------------------------------
% Test Code IGNORE
% -------------------------------------------------------------------------
% Matlab Moving Average?
% windowSize = 17;
% b = (1/windowSize)*ones(1,windowSize);
% a = 1;
% y = filter(b, a, X_mag_avg);