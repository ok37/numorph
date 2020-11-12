%% Template to run Tissue Clearing Processing Pipeline
% These are the key parameters
% Set flags to indicate how/whether to run process
adjust_intensity = ["true","load"];               % true, load, false; Whether to calculate and apply any of the following intensity adjustments. Intensity adjustment measurements should typically be performed on raw images
adjust_tile_shading = "manual";          % basic, manual, false; Can be 1xn_channels. Perform shading correction using BaSIC algorithm or using manual measurements from UMII microscope
adjust_tile_position = "true";           % true, false; Can be 1xn_channels. Normalize tile intensities by position using overlapping regions

channel_alignment = "translation";          % elastix, translation, false; Channel alignment by rigid, 2D translation or non-rigid B-splines using elastix

stitch_img = "false";                    % true, load, false; 2D iterative stitching

use_processed_images = "false";         % false or name of sub-directory in output directory (i.e. aligned, stitched...); Direct pipeline to load previously processed images in output directory
save_samples = "true";                  % true, false; Save sample results for each major step

%% Intensity Adjustment Parameters
darkfield_intensity = 101;              % 1xn_channels; Constant darkfield intensity value (i.e. average intensity of image with nothing present)
update_intensity_channels = [];         % integers; Update intensity adjustments only to certain channels 

% Parameters for manual shading corretion (i.e. adjust_tile_shading = "manual")
single_sheet = "true";                  % true, false; Whether a single sheet was used for acquisition
ls_width = [50 50 50];                  % 1xn_channels interger; Light sheet width setting for UltraMicroscope II as percentage
laser_y_displacement = 0;               % [-0.5,0.5]; Displacement of light-sheet along y axis. Value of 0.5 means light-sheet center is positioned at the top of the image

% Parameters for BaSIC shading_correction (i.e. adjust_tile_shading = "basic") 
sampling_frequency = 0.2;               % [0,1]; Fraction of images to read and sample from. Setting to 1 means use all images
shading_correction_tiles = 1:4;         % integer vector; Subset tile positions for calculating shading correction (row major order). It's recommended that bright regions are avoid
shading_smoothness = 5;                 % numeric; Factor for adjusting smoothness of shading correction. Set >1 for smoother flatfield image

%% Z Alignment Parameters
% Used for stitching and alignment by translation steps
update_z_adjustment = "false";           % true, false; Update z adjusment steps with new parameters. Otherwise pipeline will search for previously calculated parameters
z_positions = 2;                        % integer; Sampling positions along adjacent image stacks to determine z displacement. Set to 0 for no adjustment, only if you're confident tiles are aligned along z dimension
z_window = 3;                           % integer; Search window for finding corresponding tiles (i.e. +/-n z positions)
z_initial = [0 0 0];                    % 1xn_channels-1 interger; Predicted initial z displacement between reference channel and secondary channel (i.e. 

%% Channel Alignment Parameters
load_alignment_params = "true";          % true, update, false; True: apply previously calculated parameters to align individual tiles during stitching. Update: update previosuly calculated alignment parameters for specified images based on new settings
align_tiles = [];                        % Option to align only certain stacks and not all stacks. Row-major order
align_channels = [];                     % Option to align only certain channels (set to >1)
align_slices = {};                       % Option to align only certain slice ranges. Set as cell array for non-continuous ranges (i.e. {1:100,200:300})

% Specific to translation method
align_stepsize = 10;                     % interger; Only for alignment by translation. Number of images sampled for determining translations. Images in between are interpolated
save_aligned_images = "true";           % true or false. Save aligned images. If using elastix or aligning different resolutions, images will automatically be saved

% Specific to elastix method
align_chunks = [];                      % Only for alignment by elastix. Option to align only certain chunks
max_chunk_size = 300;                   % integer; Chunk size for elastix alignment. Decreasing may improve precision but can give spurious results
chunk_pad = 30;                         % integer; Padding around chunks. Should be set to value greater than the maximum expected translation in z
param_folder = "32_bins";               % 1xn_channels-1 string; Name of folders containing elastix registration parameters. Place in /supplementary_data/elastix_parameter_files/channel_alignment
mask_int_threshold = [];                % numeric; Mask intensity threshold for choosing signal pixels in elastix channel alignment. Leave empty to calculate automatically
resample_s = [3 3 1];                   % 1x3 integer. Amount of downsampling along each axis. Some downsampling, ideally close to isotropic resolution, is recommended
hist_match = [];                        % 1xn_channels-1 interger; Match histogram bins to reference channel? If so, specify number of bins. Otherwise leave empty or set to 0. This can be useful for low contrast images

%% Stitching Parameters
% Parameters for running iterative 2D stiching
stitch_sub_stack = [];                              % z positions; If only stitching a cetrain z range from all the images
stitch_sub_channel = [];                            % channel index; If only stitching certain channels
overlap = 0.15;                                     % [0,1]; overlap between tiles as fraction
blending_method = "sigmoid";                        % sigmoid, max. Recommended: sigmoid. Tile blending method
sd = 3;                                             % numeric >= 1; Steepness of sigmoid-based blending (closer to 1 gives more linear while high values (>50) is closer to block-face)
border_pad = 50;                                    % numeric >= 0; Crops borders during stitching. Increase if images shift significantly during alignment to prevent zeros values from entering stitched image
sift_refinement = "false";                          % true, false; Refine stitching using SIFT algorithm (requires vl_fleat toolbox)
use_middle = "true";                               % true, false; Recommended: false. Start stitching from the slice or optimize based on image features

%% Additional Filters and Adjustments That May Be Useful But Are Not Required
% These all currently occur during stitching step after the completed image
% has been merged.
% Parameters for rescale_intensities
rescale_intensities = "false";           % true, false; Rescaling intensities and applying gamma
lowerThresh = [100 100];                % 1xn_channels numeric; Intensity values of dimmest feature of interest. If left empty, pipeline will automatically calculate this for each channel as it is required for elastix alignment and stitching
upperThresh = [1000 2000];              % 1xn_channels numeric; Max intensity of brightest features.If left empty, pipeline will automatically calculate this for each channel as it is required for elastix channel alignment
gamma = [];                             % 1xn_channels numeric; Gamma intensity adjustment

% Background subtraction
subtract_background = "false";          % true, false. Subtrat background (similar to Fiji's rolling ball background subtraction)
nuc_radius = 13;                        % numeric >= 1; Max radius of cell nuclei along x/y in pixels. Required also for DoG filtering

% Smoothing filters
smooth_img = "false";                   % 1xn_channels, "gaussian", "median", "guided". Apply a smoothing filter
smooth_sigma = [];                      % 1xn_channels numeric; Size of smoothing kernel. For median and guided filters, it is the dimension of the kernel size

% Difference-of-Gaussian Filter
DoG_img = "false";                          % true,false; Apply difference of gaussian enhancement of blobs
DoG_minmax = [0.8,2];                       % 1x2 numeric; Min/max sigma values to take differene from.
DoG_factor = 1;                             % [0,1]; Factor controlling amount of adjustment to apply. Set to 1 for absolute DoG

% Permute image
flip_axis = "none";               % "horizontal", "vertical", "both"; Flip image along horizontal or vertical axis
rotate_axis = 0;                        % 90 or -90; Rotate image
