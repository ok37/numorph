function [lowerThresh, upperThresh, signalThresh, t_adj, y_adj, flatfield, darkfield] = measure_images(config, stack, channel_idx, get_thresholds)
%--------------------------------------------------------------------------
% Calculate various intensity adjustments including tile differences,
% light-sheet width correction, flatfield + darkfield shading correction
%--------------------------------------------------------------------------
% Some additional defaults for tile position adjustment
defaults.low_prct = 5;   % Low percentile for sampling background pixels
defaults.high_prct = 99.5; % High percentile for sampling bright pixels
defaults.pads = 0.1;     % Crop this fraction of image from the sides of images. Note this should be <0.5
defaults.image_sampling = 0.01;   % Fraction of all images to sample

% Just for measuring thresholds
if nargin <4
    get_thresholds = false;
end

% Load config parameters
resolution = config.resolution{channel_idx};

% Make sure only selected channel is chosen
stack = stack(stack.markers == config.markers(channel_idx),:);
if isempty(stack)
    error("No images found in path table for marker %s",config.markers(channel_idx));
end

% Count number of images and measure image dimensions
x_tiles = length(unique(stack.x));
y_tiles = length(unique(stack.y));

% Read image size
tempI = imread(stack.file{1});
[nrows, ncols] = size(tempI);

% Get intensity thresholds
[lowerThresh, upperThresh, signalThresh] = measure_thresholds(stack,defaults);

% If just getting thresholds, measure those and return
if get_thresholds
    return
end

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
    [flatfield, darkfield] = estimate_flatfield(config, stack);
    flatfield = single(flatfield);
    darkfield = single(darkfield);
else
    flatfield = ones(nrows,ncols,'single'); darkfield = ones(nrows,ncols,'single');
end

% Tile intensity adjustments
if x_tiles*y_tiles>1 && isequal(config.adjust_tile_position(channel_idx),"true")
    % Measure overlapping tiles and get thresholds
    defaults.image_sampling = defaults.image_sampling*10;
    [adj_matrix1,adj_matrix2] = adjust_tile_multi(config,stack,defaults);
    
    % Store tile adjustment
    t_adj(:,:,1) = adj_matrix1;
    t_adj(:,:,2) = adj_matrix2;
else
    if x_tiles*y_tiles == 1
        fprintf('%s\t Only 1 tile detected for marker %d \n',...
            datetime('now'),config.markers(channel_idx));
    end
    % Store default tile adjustments
    t_adj = [];
end

end
