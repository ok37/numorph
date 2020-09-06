%% Template to run Tissue Clearing Processing Pipeline
% Set flags to indicate how/whether to run process
use_processed_images = "false";         % Use images from aligned 		

adjust_intensity = "load";               % true, update, load, false. Intensity adjustment measurements should be performed on raw images!
adjust_ls_width = "false";               % true, false. Adjust for light-sheet width
adjust_tile_intensity = "true";          % true, false. Adjust between tile differences
shading_correction = "true";             % true, false. Perform shading correction using BaSIC algorithm
rescale_intensities = "false";           % true, false. Rescaling intensities and applying gamma 

channel_alignment = "elastix";          % elastix, translation, false
load_alignment_params = "update";         % true, update, false. True: apply previously caculated parameters. Update: update alignment for certain regions

stitch_img = "load";                    % true, load, false

use_parallel = "true";                   % true, false
number_of_cores = 14;                    % set max to either the number of cores available or the number of tiles. Set to 1 for no parallel processing

%% Channel Alignment Parameters
align_stepsize = 10;                    % Only for alignment by translation. Number of images sampled for determining translations. Images in between are interpolated
z_initial = [0 0 0];                    % Predicted initial z displacement between channels (lasers)

align_tiles = [8];                       % Option to align only certain stacks and not all stacks. Row-major order
align_channels = [3];
align_chunks = [1];
align_slices = {};
h_bins = [32 16];                       % Elastix alignment histogram bins
mask_int_threshold = 0.06;              % Mask intensity threshold for choosing signal pixels in elastix channel alignment
resample_s = [3 3 1];                   % Amount of downsampling for elastix channel alignment 

save_aligned_images = "true";           % true or false. If using elastix, images must be saved prior to stitching

%% z-Alignment Parameters for Stitching
adjust_z = "matrix";                     % true, file, matrix, false
z_positions = 15;                       % Sampling positions along adjacent image stacks to determine z displacement
z_window = 5;                           % Search window for finding corresponding tiles (i.e. +/-n z positions)

%% Stitching Parameters
stitch_sub_stack = [];               % If only stitching a cetrain z range from all the images
stitch_sub_channel = [];                % If only stitching certain channels
overlap = 0.15;                         % Overlap % between tiles as decimal
sd = 3;                                 % Steepness of sigmoid-based blending (closer to 1 gives more linear while high values (>50) is closer to block-face)
border_pad = 50;                        % Crops borders during stitching. Increase if images shift significantly during alignment to prevent zeros values from entering stitched image
blending_method = ["sigmoid","sigmoid","sigmoid"];  % sigmoid, max
sift_refinement = "true";                           % true, false
use_middle = "true";                                % true, false

%% Intensity Adjustment Parameters
single_sheet = "true";                  % true, false. Whether a single sheet was used for acquisition
shading_correction_tiles = [1:9];       % Subset tile positions for calculating shading correction (row major order)
lowerThresh = [];                       % Intensity of dimmest feature of interest. Required for alignment and stitching. If left empty, pipeline will attempt to predict this for each channel
upperThresh = [];                       % Max intensity of brightest features (ref channel first). If left empty, pipeline will attempt to predict this for each channel
gamma = [1 1 1];                        % Gamma intensity adjustment

%% Additional Filters and Adjustments During Stitching
nuc_radius = 13;                                       	% Max radius of cell nuclei in pixels
subtract_img_background = "false";                      % true, false
subtract_background = ["false", "false", "false"];     	% true, false
filter = ["false", "false", "false"];                   % Apply guided filter. Check preprocessing/apply_diffuse.m for parameters
dog = [];                                               % Apply difference of gaussian enhancement of blobs
