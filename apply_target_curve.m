% DSTII Final Project %
% Harrison Zafrin hzz200
% x_t = time domain signal
% fs = sampling rate
% win_size = RMS window size, this equals FFT size in other analysis
% -------------------------------------------------------------------------
% Perform spectral matching via custom IIR filter with coefficients
% generated via least squares Yule Walker.  A spectral analysis is
% performed, and the desired magnitude resposne is obtained to create the
% desired transfer funnction
% -------------------------------------------------------------------------
function [ filtered_output ] = apply_target_curve( x_t, T_mag, fftparams,...
                                                    fs, active_frames)

% -------------------------------------------------------------------------
% Perform an STFT of the incoming signal
% -------------------------------------------------------------------------
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

% Thresholding of magnitudes below 0.0001
% X_mag = X_mag(X_mag > 0.0001);

[ X_mag ] = normalize_magMatrix( X_mag );

% -------------------------------------------------------------------------
% Obtain Desired Magnitude Response Per Frame
% -------------------------------------------------------------------------

% Array of the freqs we want, 1/3 octave center freqs, 33 freq bands
hertz_vector = [0, 20, 25, 31.5, 40, 50, 63, 80, 100, 125, 160, 200, 250,...
                315, 400, 500, 630, 800, 1000, 1250, 1600, 2000, 2500, 3150,...
                4000, 5000, 6300, 8000, 10000, 12220, 16000, fs/2];

% The bins in bin_vector correspond to the freqs in hertz_vector
bin_vector = 1+floor((fftparams.win_size-1)*(hertz_vector/fs));

% Obtain Desired Magnitude Response for a Frame
for i=1:length(X_mag)
    
    % Grab a frame
    frame = X_mag(:,i);
    
    % implement 3.3.1 dry/wet HERE!!!!!!    

    Hd_Mag(:,i) = T_mag(bin_vector)./frame(bin_vector)';

end

[ Hd_Mag ] = normalize_magMatrix( Hd_Mag );

% -------------------------------------------------------------------------
% Filter Curve Smoothing via Exponential Moving Average
% EMA = Hm'(w) = alpha* Hm-1'(w) + (1-alpha) * Hm(w)
% Hm is the current frame
% -------------------------------------------------------------------------

% Determine degree of filtering
alpha = exp^(-1/(0.5*fs));

% Create Zeroed frame for intiial m-1 in equation
zeroframe = zeros(size(Hd_Mag, 2), 1);

% Concat the two together like a zero-pad of the matrix
Hd_Mag = [zeroframe Hd_Mag];

% Run through the EMA
for i=2:length(Hd_mag)
   
    Hd_Mag(:,i) = alpha * Hd_Mag(:,i-1) + (1-alpha) * Hd_Mag(:,i);
    
end

% Remove first frame
Hd_Mag = Hd_Mag(:,2:end);

% -------------------------------------------------------------------------
% Obtain filter coefficients via Yule-Walker, N is 16 accroding to lit.
% -------------------------------------------------------------------------

[b,a] = yulewalk(16, hertz_vector/(fs/2) , Hd_Mag(:,i));

end
