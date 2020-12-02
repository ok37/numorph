function [t_adj, lowerThresh, upperThresh, y_adj, flatfield, darkfield] = measure_images(config, stack, channel_idx)
%--------------------------------------------------------------------------
% Calculate various intensity adjustments including tile differences,
% light-sheet width correction, flatfield + darkfield shading correction
%--------------------------------------------------------------------------
% Some additional defaults for tile position adjustment
defaults.low_prct = 5;   % Low percentile for sampling background pixels
defaults.high_prct = 99; % High percentile for sampling bright pixels
defaults.pads = 0.25;     % Crop this fraction of image from the sides of images
defaults.image_sampling = 0.1;   % Fraction of all images to sample

% Load config parameters
resolution = config.resolution{channel_idx};

% Count number of images and measure image dimensions
x_tiles = length(unique(stack.x));
y_tiles = length(unique(stack.y));

% Read image size
tempI = imread(stack.file{1});
[nrows, ncols] = size(tempI);

% Calculate and adjust for laser width
if isequal(config.adjust_tile_shading,"manual")
    fprintf('%s\t Adjusting For Laser Width \n',datetime('now'));    
    y_adj = adjust_ls_width_measured(config.single_sheet(channel_idx),tempI,...
    config.ls_width(channel_idx),resolution,config.laser_y_displacement(channel_idx))';
else
    y_adj = ones(nrows,1);
end  

% Calculate shading correction using BaSiC
if isequal(config.adjust_tile_shading(channel_idx),"basic")
    %[flatfield, darkfield] = estimate_flatfield(config, stack);
    flatfield = single(flatfield);
    darkfield = single(darkfield);
else
    flatfield = ones(nrows,ncols,'single'); darkfield = ones(nrows,ncols,'single');
end

% Tile intensity adjustments
if x_tiles*y_tiles>1 && isequal(config.adjust_tile_position(channel_idx),"true")
    % Measure overlapping tiles and get thresholds
    [adj_matrix1,adj_matrix2,lowerThresh,upperThresh] = adjust_tile_multi(config,stack,defaults);
    
    % Store tile adjustment
    t_adj(:,:,1) = adj_matrix1;
    t_adj(:,:,2) = adj_matrix2;
else
    if x_tiles*y_tiles == 1
        fprintf('%s\t Only 1 tile detected for channel %d \n',datetime('now'),channel_idx);
    end
    % Get thresholds without measuring overlapping tiles
    [lowerThresh,upperThresh] = adjust_tile_single(stack,defaults);
    
    % Store default tile adjustments
    t_adj = [];
end

end
