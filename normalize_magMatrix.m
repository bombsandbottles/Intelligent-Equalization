% DSTII Final Project %
% Harrison Zafrin hzz200
% x = is a framed magnitude matrix | rows are bins and columns are frames
% -------------------------------------------------------------------------
% Normalizes a magnitude matrix to (0,1) as in the literature
% -------------------------------------------------------------------------
function [ normalized_mags ] = normalize_magMatrix( x )

% Normalize Each Frame (0,1).  Divide mags by max mag per frame
for i=1:size(x, 2);
    
    % Grab a frame
    frame = x(:,i);
    
    % Block NaNs
    if max(frame) == 0
        x(:,i) = 0;
    else
        % Normalize a frame
        x(:,i) = frame/(max(frame));
    end
    
end

% LOAD OUT
normalized_mags = x;

end

