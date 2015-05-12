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
function [ filtered_output, x_t_filt, x_t_windowed ] = apply_target_curve( x_t, T_mag, fftparams,...
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

% -------------------------------------------------------------------------
% Threshold magnitudes below 0.0001
% -------------------------------------------------------------------------
for i=1:size(X_mag, 2)
    
    % Grab a frame
    frame = X_mag(:,i);
    
    for j=1:length(frame)
        if frame(j) < 0.0001
            frame(j) = 0.0001;
        end
    end
    
    % Put back into Magitude matrix    
    X_mag(:,i) = frame;
    
end

% [ X_mag ] = normalize_magMatrix( X_mag );

% -------------------------------------------------------------------------
% Obtain Desired Magnitude Response Per Frame
% -------------------------------------------------------------------------

% Array of the freqs we want, 1/3 octave center freqs, 33 freq bands
hertz_vector = [0, 20, 25, 31.5, 40, 50, 63, 80, 100, 125, 160, 200, 250,...
                315, 400, 500, 630, 800, 1000, 1250, 1600, 2000, 2500, 3150,...
                4000, 5000, 6300, 8000, 10000, 12220, 16000, fs/2];

% The bins in bin_vector correspond to the freqs in hertz_vector
bin_vector = 1+floor((fftparams.win_size-1)*(hertz_vector/fs));

T_mag = T_mag(bin_vector);

% Obtain Desired Magnitude Response
for i=1:size(X_mag, 2)
    
    % Grab a frame
    frame = X_mag(:,i);
    frame = frame';
    frame = frame(bin_vector);
    
    % Perform |T(w)| / |X(w)| for each w (frequency bin)    
    for j=1:length(frame)
        desired_mag(j) = T_mag(j) / frame(j);
    end

    % Create the Desired Magnitude Response for Frame (i)
    Hd_Mag(:,i) = desired_mag;

end

% Create previous frame
% prev_frame = zeros(size(X_mag, 1), 1);

% % Obtain Desired Magnitude Response for a Frame
% for i=1:size(X_mag, 2)
%     
%     % If the Frame is Active, We Grab it
%     if active_frames(i) == 1
%         
%         % Grab a frame
%         frame = X_mag(:,i);
%         
%     % Otherwise we use the previous active frame
%     elseif active_frames(i) == 0
%         
%         frame = prev_frame;
%         
%     end
%   
%     % implement 3.3.1 dry/wet HERE!!!!!!    
% 
%     % Create the Desired Magnitude Response for Frame (i)
%     Hd_Mag(:,i) = T_mag(bin_vector)./frame(bin_vector)';
%     
%     % Store current frame into previous frame
%     prev_frame = frame;
% 
% end

% Remove NaNs and -Infs
% Hd_Mag(isnan(Hd_Mag)) = 0;
% Hd_Mag(isinf(Hd_Mag)) = 0;


% -------------------------------------------------------------------------
% Normalize Hd_Mag so that |Hd(w)| is between 0 and 1
% -------------------------------------------------------------------------
% Noramlize Across COLUMNS, Coefficients Across Time
for i=1:size(Hd_Mag, 2)
        Hd_Mag(:, i) = ( Hd_Mag(:, i) )/ norm(Hd_Mag(:, i), 2);
end

% Hd_Mag = abs(Hd_Mag);


% -------------------------------------------------------------------------
% Filter Curve Smoothing via Exponential Moving Average
% EMA = Hm'(w) = alpha* Hm-1'(w) + (1-alpha) * Hm(w)
% Hm is the current frame
% -------------------------------------------------------------------------

% Determine degree of filtering
% alpha = exp^(-1/(0.5*fs));
alpha = 0.9;

% Make output for prev_value storage
output = 0;

% Run through the EMA
for i=2:size(Hd_Mag, 2)
   
    Hd_Mag(:,i) = (alpha * output) + ((1-alpha) * Hd_Mag(:,i));
    output = Hd_Mag(:,i);
end

% -------------------------------------------------------------------------
% Obtain filter coefficients via Yule-Walker, N is 16 accroding to lit.
% Then filter with an IIR filter
% -------------------------------------------------------------------------

% Create Normalized Frequency Vector
f = hertz_vector/(fs/2);
f = f';

% Pre-Allocate Zeros for Difference Equation
y_buff = zeros(1,16); 
x_buff = zeros(1,16);

% Filter x_t frame by frame with dynamic magnitude via Yule-Walker
for i=1:size(Hd_Mag, 2)
    
    % Get IIR coeffs for a frame
    [b, a] = yulewalk(16, f, Hd_Mag(:,i));
    
    % Grab active frame
    frame = x_t_windowed(:,i);
    
    % Apply IIR filter via difference equation
    for n=1:length(frame)
       
        % Run the diff equation
        filt_frame(n) = (b(1)*frame(n))     + (b(2)*x_buff(1))   + ...
                        (b(3)*x_buff(2))    + (b(4)*x_buff(3))   + ...
                        (b(5)*x_buff(4))    + (b(6)*x_buff(5))   + ...
                        (b(7)*x_buff(6))    + (b(8)*x_buff(7))   + ...
                        (b(9)*x_buff(8))    + (b(10)*x_buff(9))  + ...
                        (b(11)*x_buff(10))  + (b(12)*x_buff(11)) + ...
                        (b(13)*x_buff(12))  + (b(14)*x_buff(13)) + ...
                        (b(15)*x_buff(14))  + (b(16)*x_buff(15)) + ...
                        (b(17)*x_buff(16))  - ...
                        (a(2)*y_buff(1))    - ...
                        (a(3)*y_buff(2))    - (a(4)*y_buff(3))   - ...
                        (a(5)*y_buff(4))    - (a(6)*y_buff(5))   - ...
                        (a(7)*y_buff(6))    - (a(8)*y_buff(7))   - ...
                        (a(9)*y_buff(8))    - (a(10)*y_buff(9))  - ...
                        (a(11)*y_buff(10))  - (a(12)*y_buff(11)) - ...
                        (a(13)*y_buff(12))  - (a(14)*y_buff(13)) - ...
                        (a(15)*y_buff(14))  - (a(16)*y_buff(15)) - ...
                        (a(17)*y_buff(16));
        
        % Shift the Samples in the Equation So That x-1 == n etc.
        x_buff(2:end) = x_buff(1:end-1);
        x_buff(1) = frame(n);
        
        % Shift the Samples in the Equation So That y-1 == n etc.
        y_buff(2:end) = y_buff(1:end-1);
        y_buff(1) = filt_frame(n);
        
    end
    
    % Load into Filtered Matrix
    x_t_filt(:,i) = filt_frame;
    
end

% -------------------------------------------------------------------------
% OLA Unbuffer, This section converts the filtered matrix back to mono 
% audio stream
% -------------------------------------------------------------------------

% Pre-Allocate Output Vector
[win_length, num_win] = size(x_t_filt);
hop_size = win_length/2;
x_t_unbuffered = zeros(((num_win-1) * hop_size)+win_length, 1);

% Unbuffer the Windowed Kernel
for i=1:num_win
    seg_start    = (i-1)*hop_size+1;
    seg_end      = (i*win_length)-(hop_size*(i-1));
    x_t_unbuffered(seg_start:seg_end)  = x_t_unbuffered(seg_start:seg_end) + x_t_filt(:,i); 
end

% Transpose the Signal
filtered_output = x_t_unbuffered';

end
