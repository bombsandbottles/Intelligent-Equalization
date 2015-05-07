% DSTII Final Project %
% Harrison Zafrin hzz200
% x_t = time domain signal
% fs = sampling rate
% win_size = RMS window size, this equals FFT size in other analysis
% -------------------------------------------------------------------------
% Loudness Calculation for EBU loudness R- 128
% -------------------------------------------------------------------------
function [ LU, active_frames ] = calc_loudness_EBU( x_t, fs, fftparams )

% Resample for K-Filter, coefficients are made for 48khz
% [p, q] = rat(48000/fs);
% x_t = resample(x_t, p, q);

% -------------------------------------------------------------------------
% K-frequency weighting (perceptual filtering)
% -------------------------------------------------------------------------

% Create Stage 1 of Pre-Filter, High Shelf Boost of +4dB at 1681hz
a1=[1,-1.69065929318241,0.73248077421585];
b1=[1.53512485958697,-2.69169618940638,1.19839281085285];
    
x_t_hishelf = filter(b1, a1, x_t);

% Create Stage 2 of Pre-Filter, high pass at 38hz
a2=[1,-1.99004745483398,0.99007225036621];
b2=[1.0,-2.0,1.0];
    
x_t_kfilt = filter(b2, a2 ,x_t_hishelf);

% -------------------------------------------------------------------------
% Perform Modified RMS to get LU measurement per window
% -------------------------------------------------------------------------

% Get amount of sample overlap per window
n_overlap = fftparams.win_size - fftparams.hop_size;

% Create Window a hanning window to prevent spectral leakage
window = hann(fftparams.win_size);

% Buffer with 4096 overlap to match FFT frames
x_t_buff = buffer(x_t_kfilt, fftparams.win_size, n_overlap, 'nodelay');

% Create Window Matrix
window_mat = repmat(window, 1, size(x_t_buff, 2));

% Window the Signal
x_t_windowed = x_t_buff .* window_mat;

% Square each element, mean over the window size, convert to dB scale
LU = 0.691 * 10*log10(mean(x_t_windowed.^2)+eps);

% Apply an EMA filter to provide estimation of energy over 3-seconds
% Y[n] = alpha * Y[n-1] + (1-alpha) * X[n]
alpha = 0.9;
output = 0;
for i=1:length(LU)
    LU(i) = (alpha*output) + ((1-alpha) * LU(i));
    output = LU(i);
end


% -------------------------------------------------------------------------
% Apply Hysteresis Noise Gate to Signal to Determine Active Frames
% -------------------------------------------------------------------------

% Define Upper and Lower Thresholds in -25LUFS and -30LUFS
upper_thresh = -25;
lower_thresh = -30;

% Pre-Allocate Active Frames
active_frames = zeros(1, length(LU));

% Past Frame = 0 At first
past_frame = 0;

for i=1:length(LU)
    
    % If above -25LUFS, 1 equal active frame
    if LU(i) > upper_thresh
        active_frames(i) = 1;

    % If the past frame was above -25LUFS and this ones above -30
    elseif past_frame > upper_thresh && LU(i) > lower_thresh
        active_frames(i) = 1;
        
    % If the current frame is below the threshold, non-active
    elseif LU(i) < upper_thresh
        active_frames(i) = 0;
        
    end
    
    % The analyzed frame becomes the past frame
    past_frame = LU(i);
    
end

end