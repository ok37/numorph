%% Template to run Tissue Clearing Processing Pipeline
% These are the key parameters
% Set flags to indicate how/whether to run process
adjust_intensity = "false";             % true, update, load, false; Whether to calculate and apply intensity adjustments. Intensity adjustment measurements should be performed on raw images!
shading_correction = "true";            % true, false; Perform shading correction using BaSIC algorithm
adjust_ls_width = "false";              % true, false; Adjust for light-sheet width using manual measurements
adjust_tile_intensity = "true";         % true, false; Adjust between tile differences

channel_alignment = "elastix";          % elastix, translation, false; Channel alignment by rigid, 2D translation or non-rigid B-splines using elastix

stitch_img = "load";                    % true, load, false; 2D iterative stitching

use_processed_images = "false";         % false or name of sub-directory in output directory (i.e. aligned, stitched...); Direct pipeline to load previously processed images in output directory
save_samples = "true";                  % true, false; Save sample results for each major step

%% Intensity Adjustment Parameters
% Parameters for adjust_ls_width
single_sheet = "true";                  % true, false; Whether a single sheet was used for acquisition
ls_width = [50 50 50];                  % 1xn_channels interger value; Light sheet width setting for UMII as percentage
laser_y_displacement = 0;               % [-0.5,0.5]; Rough displacement of light-sheet along y axis. Value of 0.5 means light-sheet center is positioned at the top of the image

% Parameters for shading_correction
sampling_frequency = 0.2;               % [0,1]; Fraction of images to read and sample from 
shading_correction_tiles = 1:9;         % integer vector; Subset tile positions for calculating shading correction (row major order)
shading_smoothness = 1.5;               % numeric; Factor for adjusting smoothness of shading correction. Set >1 for smoother flatfield image

%% Z Alignment Parameters
% Used for stitching and alignment by translation steps
update_z_adjustment = "false";          % true, false; Update z adjusment steps with new parameters. Otherwise pipeline will search for previously calculated parameters
z_positions = 2;                        % interger value; Sampling positions along adjacent image stacks to determine z displacement. Set to 0 for no adjustment, only if you're confident tiles are aligned along z dimension
z_window = 5;                           % interger value; Search window for finding corresponding tiles (i.e. +/-n z positions)
z_initial = [0 0 0];                    % 1xn_channels interger value; Predicted initial z displacement between channels (lasers)

%% Channel Alignment Parameters
load_alignment_params = "true";          % true, update, false; True: apply previously caculated parameters to align individual tiles. Update: update previosuly calculated alignment parameters for specified images based on new settings
align_tiles = [];                        % Option to align only certain stacks and not all stacks. Row-major order
align_channels = [];                     % Option to align only certain channels (set to >1)
align_slices = {};                       % Option to align only certain slice ranges. Set as cell array for non-continuous ranges (i.e. {1:100,200:300})
save_aligned_images = "false";           % true or false. Save aligned images. If using elastix, images must be saved prior to stitching

% Specific to translation method
align_stepsize = 10;                    % Only for alignment by translation. Number of images sampled for determining translations. Images in between are interpolated

% Specific to elastix method
align_chunks = [];                      % Only for alignment by elastix. Option to align only certain chunks
max_chunk_size = 300;                   % Chunk size for elastix alignment. Decreasing may improve precision but can give spurious results
chunk_pad = 30;                         % Padding around chunks. Should be set to value greater than the maximum expected translation in z
h_bins = [32 16];                       % Elastix alignment histogram bins. Set for n-1 markers to be registered
mask_int_threshold = [];                % Mask intensity threshold for choosing signal pixels in elastix channel alignment. Leave empty to calculate automatically
resample_s = [3 3 1];                   % 1x3 interger. Amount of downsampling along each axis. Some downsampling, ideally close to isotropic resolution, is recommended
hist_match = [0 0 128];             % Match histogram bins to reference channel? If so, specify number of bins. Otherwise leave at 0. This can be useful for low contrast images

%% Stitching Parameters
% Parameters for running iterative 2D stiching
stitch_sub_stack = [];                              % z positions; If only stitching a cetrain z range from all the images
stitch_sub_channel = [];                            % channel index; If only stitching certain channels
overlap = 0.15;                                     % [0,1]; overlap between tiles as fraction
blending_method = "sigmoid";                        % sigmoid, max. Recommended: sigmoid. Tile blending method
sd = 3;                                             % numeric >= 1; Steepness of sigmoid-based blending (closer to 1 gives more linear while high values (>50) is closer to block-face)
border_pad = 50;                                    % numeric >= 0; Crops borders during stitching. Increase if images shift significantly during alignment to prevent zeros values from entering stitched image
sift_refinement = "false";                          % true, false; Refine stitching using SIFT algorithm (requires vl_fleat toolbox)
use_middle = "false";                               % true, false; Recommended: false. Start stitching from the slice or optimize based on image features

%% Additional Filters and Adjustments That May Be Useful But Are Not Required
% Parameters for rescale_intensities
rescale_intensities = "false";          % true, false; Rescaling intensities and applying gamma
lowerThresh = [];                       % 1xn_channels numeric; Intensity values of dimmest feature of interest. Required for elastix alignment and stitching. If left empty, pipeline will automatically calculate this for each channel
upperThresh = [];                       % 1xn_channels numeric; Max intensity of brightest features. Required for elastix channel alignment. If left empty, pipeline will automatically calculate this for each channel
gamma = [1 1 1];                        % 1xn_channels numeric; Gamma intensity adjustment

nuc_radius = 13;                                       	% numeric >= 1; Max radius of cell nuclei in pixels
subtract_background = ["false", "false", "false"];     	% true, false. Subtrat background (similar to Fiji's rolling ball background subtraction)
filter = ["false", "false", "false"];                   % Apply guided filter. Check preprocessing/apply_diffuse.m for parameters
dog = [];                                               % Apply difference of gaussian enhancement of blobs
