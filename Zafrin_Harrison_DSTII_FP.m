% DSTII Final Project %
% Harrison Zafrin hzz200
% Intelligent Equalization tool using Yule-Walker Method
% -------------------------------------------------------------------------
clear;
clc;

% -------------------------------------------------------------------------
% Create Target Database
% -------------------------------------------------------------------------
% Female Pop Vocal References
filenames = {'BS_BTI.wav',...
             'CA_YB.wav',...
             'KP_FW.wav'...
             'KP_ET.wav'...
             'K_DY.wav'};

% -------------------------------------------------------------------------
% Create fft params structure to pass to all functions 
% -------------------------------------------------------------------------
% FFT Params for use in Loudness Calculation
field1 = 'win_size'; win_size = 4096;
field2 = 'hop_size'; hop_size = win_size/2;

% Create ffrparams Structure
fftparams = struct(field1, win_size,...
                   field2, hop_size);

% -------------------------------------------------------------------------
% Compute Target Equalization Curve
% -------------------------------------------------------------------------

% Compute the magnitude spectra for each song in the dataset
for i=1:length(filenames)

    % Obtain normalized, averaged 4096 point magnitude spectra for a file
    [ X_mag, X_mag_mean, X_mag_cum, fs ] = average_spectra( filenames{i}, fftparams );
    
    % Load the mean magnitudes and cumululative normalized magnitudes
    mag_mean_matrix(:,i) = X_mag_mean;
    mag_cum_matrix(i,:) = X_mag_cum;
    
end

% Compute the average magnitude spectra across all songs in the dataset
[ X_mag_avg, freq_vector ] = create_overall_average_spectrum( mag_mean_matrix,...
                                                              mag_cum_matrix,...
                                                              size(mag_cum_matrix, 1),...
                                                              fs);
% Convert to dB
X_mag_avg = mag2db(X_mag_avg);

% Smooth the curve with a 17 point moving average filter using MATLAB filt
[ T_mag_matlab ] = movingavg( X_mag_avg, 17, fftparams, fs );

% Smooth the curve with a 17 point moving average filter
[ T_mag ] = movingavgfilter_17pnt( X_mag_avg, 17, fftparams, fs );

% -------------------------------------------------------------------------
% Figure comparing averaged Target EQ Curves
% -------------------------------------------------------------------------
figure;
semilogx(X_mag_avg);
hold on;
semilogx(T_mag_matlab, 'r');
hold on;
semilogx(T_mag, 'g');
set(gca,'XTickLabel',num2str(get(gca,'XTick').'));

% -------------------------------------------------------------------------
% Apply Target Curve to Analyzed Audio
% -------------------------------------------------------------------------

% Import Audio File
[x_t, fs, t] = import_audio('Three Nineteen Fifteen.aif');

% Noise Test
% fs = 44100;
% x_t = rand(1,44100*10); 
% x_t  = x_t - mean(x_t);

% Determine Active Frames for Analysis
[ LU, active_frames ] = calc_loudness_EBU( x_t, fs, fftparams );

% Filter the audio
[ filtered_output, x_t_filt ] = apply_target_curve( x_t, T_mag, fftparams, fs, active_frames );
                                                          
% -------------------------------------------------------------------------
% Test Plot
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
