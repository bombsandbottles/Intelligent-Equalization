function [x_t, fs, t] = import_audio(filepath)
% Write a Matlab function that reads any audio (wav) file, and returns the signal x(t), the
% vector of time points t in seconds, and the sample rate fs. If the audio is stereo, the second
% channel should be disregarded. (This will be a very short function). [0.5 pts]

%   Read Audio File, Return fs and Signal
[x_t, fs] = audioread(filepath);

%   If the Audio File is Stereo, Discard The Right Channel and Make Mono
if size(x_t,2) > 1
	x_t = x_t(:,1);
end

%   Transpose to 1 x T array
x_t = x_t';

%   Create Time Vector t
time = length(x_t)/fs;
t = linspace(0, time, length(x_t));


end

