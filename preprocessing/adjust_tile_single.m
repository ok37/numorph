function [lowerThresh,upperThresh] = adjust_tile_single(stack,defaults)
%--------------------------------------------------------------------------
% Calculate intensity adjustment thresholds for single tile layout by 
% measuring whole images.
%--------------------------------------------------------------------------

% Load defaults
low_prct = defaults.low_prct;   % Low percentile for sampling background pixels
high_prct = defaults.high_prct; % High percentile for sampling bright pixels
image_sampling = defaults.image_sampling;   % Fraction of all images to sample

% Get vector of evenly psaced positions
stack = sortrows(stack,{'x','y','z'});
nb_images = round(height(stack)*image_sampling);
z_pos = round(linspace(1,height(stack),nb_images));

% Read images
tempI = imread(stack.file{1});
dims = size(tempI);
I = zeros([dims nb_images],'uint16');
I(:,:,1) = tempI;
for i = 2:length(z_pos)
    I(:,:,i) = imread(stack.file{z_pos(i)});
end

% Get percentile values
lowerThresh = prctile(I(:),low_prct);
upperThresh = prctile(I(:),high_prct);
fprintf('%s\t Measured Lower and Upper Intensities:\t %.1f\t %.1f\n',datetime('now'),...
    lowerThresh,upperThresh);  

end