function [lowerThresh, upperThresh] = get_thresholds(stack)
% Get upper and lower intensity thresholds from a small subset of images in
% a stack
sampling_prctile_low = 1;
sampling_prctile_high = 99;
sampling_freq = 0.01;

% Get vector of evenly psaced positions
stack = sortrows(stack,{'x','y','z'});
n_images = round(height(stack)*sampling_freq);
z_pos = round(linspace(1,height(stack),n_images));

% Read images
tempI = imread(stack.file{1});
dims = size(tempI);
I = zeros([dims n_images],'uint16');
I(:,:,1) = tempI;
for i = 2:length(z_pos)
    I(:,:,i) = imread(stack.file{z_pos(i)});
end

% Get percentile values
lowerThresh = prctile(I(:),sampling_prctile_low);
upperThresh = prctile(I(:),sampling_prctile_high);
end
    
    





