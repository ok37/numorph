%% Template to run Tissue Clearing Processing Pipeline
% Set flags to indicate how/whether to run process
use_processed_images = "false";         % false or name of sub-directory in output directory (i.e. aligned, stitched...)
save_samples = "true";                  % true, false. Save sample results for each major step

adjust_intensity = "load";              % true, update, load, false. Intensity adjustment measurements should be performed on raw images!
adjust_ls_width = "true";               % true, false. Adjust for light-sheet width using manual measurements
adjust_tile_intensity = "true";         % true, false. Adjust between tile differences
shading_correction = "false";           % true, false. Perform shading correction using BaSIC algorithm
rescale_intensities = "true";           % true, false. Recommended: false. Rescaling intensities and applying gamma.

channel_alignment = "translation";      % elastix, translation, false

stitch_img = "load";                    % true, load, false

%% Intensity Adjustment Parameters
% Parameters for adjust_ls_width
single_sheet = "true";                  % true, false. Whether a single sheet was used for acquisition
ls_width = [50 50 50];                  % light sheet width
laser_y_displacement = 0;               % displacement of light-sheet along y axis

% Parameters for shading_correction
sampling_frequency = 0.2;               % Fraction of image to read and sample from 
shading_correction_tiles = 1:9;         % Subset tile positions for calculating shading correction (row major order)
shading_smoothness = 1.5;               % Factor to adjust shading correction smoothness. Set >1 for smoother flatfield image

% Parameters for rescale_intensities
lowerThresh = [];                       % Intensity of dimmest feature of interest. Required for alignment and stitching. If left empty, pipeline will attempt to predict this for each channel
upperThresh = [];                       % Max intensity of brightest features (ref channel first). If left empty, pipeline will attempt to predict this for each channel
gamma = [1 1 1];                        % Gamma intensity adjustment

%% Z Alignment Parameters
% Used for stitching and alignment by translation steps
update_z_adjustment = "false";          % true, false. Update z adjusment steps with new parameters. Otherwise pipeline will search for previously calculated parameters
z_positions = 2;                        % Sampling positions along adjacent image stacks to determine z displacement. Set to 0 for no adjustment, only if you're confident tiles are aligned along z dimension
z_window = 5;                           % Search window for finding corresponding tiles (i.e. +/-n z positions)
z_initial = [0 0 0];                    % Predicted initial z displacement between channels (lasers)

%% Channel Alignment Parameters
load_alignment_params = "true";         % true, update, false. True: apply previously caculated parameters for individual tiles. Update: update alignment for certain regions
align_tiles = [8];                      % Option to align only certain stacks and not all stacks. Row-major order
align_channels = [3];                   % Option to align only certain channels (set to >1)
align_slices = {};                      % Option to align only certain slice ranges. Set as cell array for non-continuous ranges (i.e. {1:100,200:300})
save_aligned_images = "true";           % true or false. Save aligned images. If using elastix, images must be saved prior to stitching

% Specific to translation method
align_stepsize = 10;                    % Only for alignment by translation. Number of images sampled for determining translations. Images in between are interpolated

% Specific to elastix method
align_chunks = 1;                       % Only for alignment by elastix. Option to align only certain chunks
h_bins = [32 16];                       % Elastix alignment histogram bins. Set for n-1 markers to be registered
mask_int_threshold = [];                % Mask intensity threshold for choosing signal pixels in elastix channel alignment. Leave empty to calculate automatically
resample_s = [3 3 1];                   % Amount of downsampling for elastix channel alignment 

%% Stitching Parameters
stitch_sub_stack = [];                              % z positions; If only stitching a cetrain z range from all the images
stitch_sub_channel = [];                            % channel index; If only stitching certain channels
overlap = 0.15;                                     % overlap between tiles as decimal
blending_method = ["sigmoid","sigmoid","sigmoid"];  % sigmoid, max
sd = 3;                                             % Steepness of sigmoid-based blending (closer to 1 gives more linear while high values (>50) is closer to block-face)
border_pad = 50;                                    % Crops borders during stitching. Increase if images shift significantly during alignment to prevent zeros values from entering stitched image
sift_refinement = "true";                           % true, false. Refine stitching using SIFT algorithm (requires vl_fleat toolbox)
use_middle = "true";                                % true, false. Recommended: false. Start stitching from the slice or optimize based on image features

%% Additional Filters and Adjustments During Stitching
nuc_radius = 13;                                       	% Max radius of cell nuclei in pixels
subtract_img_background = "false";                      % true, false. Optimized background subtraction (see 
subtract_background = ["false", "false", "false"];     	% true, false
filter = ["false", "false", "false"];                   % Apply guided filter. Check preprocessing/apply_diffuse.m for parameters
dog = [];                                               % Apply difference of gaussian enhancement of blobs
