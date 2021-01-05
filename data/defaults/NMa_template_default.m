%% Template to run Tissue Clearing Analysis Pipeline
% These are the key parameters
% Set flags to indicate how/whether to run process

% Registration
register_image = "true";         % true, load, false; Register image to reference atlas. Option 'load' will load registration parameters
generate_mask = "false";         % true, load, false; Generate annotation using Allen ccf v3 annotations 

% Nuclei detection
count_cells = "3dunet";               % 3dunet, hessian, load
use_mask = "true";                    % true, false; Use mask for cell counting
measure_local_background = "load";    % true, false; measure local image background

% Cell-type classification
count_colocalized = "false";          % gmm, svm, threshold, false; Cell-type classification method

use_processed_images = "stitched";         % false or name of sub-directory in output directory (i.e. aligned, stitched...); Direct pipeline to load previously processed images in output directory
save_counts = "true";                   % true, false, overwrite 

%% Registration Parameters
% Resampling
ara_resolution = [25,25,25];            % Resolution of Allen Reference Atlas
recalculate_resampling = "false";       % true, false
mask_file = [];                         % Name of custom defined mask file. Otherwise, 
mask_resolution = [];                   % Image resolution of mask file

% Registration
registration_method = "p";                  % Affine, BSpline, Points, Other
structures_of_interest = "cortex";          % Structure of interest
atlas_file = "ara_nissl_left_10.nii";       % Name of the atlas file to register to. Can provide full path or place .nii file in supplementary_data
hemisphere = "left";                        % left, right, whole (which brain hemisphere or whole brain)
orientation = 'lateral';                    % lateral, dorsal, ventral. Which orientation is at z=0 (i.e. for hemisphere laying flat on midline, this would be lat
direction = 'img_to_atlas';                 % Forward = register image to atlas; Reverse = register atlas to image
calculate_inverse = 'true';                 % true, false. Whether to calculate the inverse transform using elastix's Displacement Magnitude Penalty
mask_coordinates = [];                      % Row start, row end, col start, col end
save_registered_image = 'true';             % Save a copy of registration results

%% Nuclei Detection
% 3D-Unet
model_file = '128_model.h5';    % Model file name located in /analysis/3dunet/nuclei/models
chunk_size = [112,112,32];      % Chunk size of unet model
chunk_overlap = [16,16,8];      % Overlap between chunks

trained_resolution = [1.21,1.21,4]; % Resolution at which the model was trained. Only required if resampling
resample_chunks = "false";          % Resample images to match model resolution. This process takes significantly longer

gpu = '0';                      % Cuda visible device index

% Hessian blob detector
min_intensity = 200;                    % Minimum intensity for cell nuclei
nuc_diameter_range = [6,20];            % Range of nuclei diameters (in pixels)

%% Cell-Type Classification
z_normalization = 'false';                  % Apply z normalization to thresholds
log_outliers = 0;                           % Apply log to z scores above or below this threshold
skip_c1 = "true";                           % Skip classifying 1st (reference/nuclei) channel

% Threshold classification
thresholds = [0.005,0.005];                 % Intensity thresholds 
expression{1} = "1*mode + 5*mad";           % Threshold expression
expression{2} = "1*mode + 2*mad";           % Threshold expression

% GMM classification
n_clusters = 4;                             % Number of clusters for GMM
mix_proportions = [0.55,0.15,0.30,0.01];	% Intial mixing proportions for n-1 markers for intializing GMM. Leave empty to estimate with k-means++
mix_markers = {[],2,3,[2,3]};               % Marker numbers for the mixing proprtions described above
confidence = [0.5,0.5,0.5,0.5];             % Posterior probability of cell being positive for a marker. This will get adjusted with addition of more clusters
stratify_structures = "true";               % Run GMM clustering on individual structures
split_markers = "false";                    % Run GMM seperately for each marker

% SVM classification


