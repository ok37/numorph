%% Template to run Tissue Clearing Analysis Pipeline
% These are the key parameters
% Set flags to indicate how/whether to run process
resample_images = "true";                       % true, update, false; Perform image resampling 
register_images = "false";                      % true, update, false; Register image to reference atlas.
count_nuclei = "false";                         % true, update, false; Count cell nuclei or other blob objects.
classify_cells = "false";                       % true, update, false; Classify cell-types for detected nuclei centroids.

use_processed_images = "stitched";              % false or name of sub-directory in output directory (i.e. aligned, stitched...); Direct pipeline to load previously processed images in output directory
structures = "structure_template.csv";          % Specify csv file in ./annotations detailing which structures to analyze
save_samples = "true";                          % true, false

%% Resampling Parameters
% Resampling
resample_channels = [];                     % Resample specific channels.
resample_resolution = 25;                   % Resample resolution (. Should match resolution of Allen Reference Atlas

%% Registration Parameters
registration_parameters = [];                % Name of folder containing elastix registration parameters. Place in /data/elastix_parameter_files/atlas_registration. If empty, will load from "default" or "points"
register_channels = 1;                       % Which channel to register to atlas
atlas_file = "ara_nissl_25.nii";             % Name of the atlas file to register to. Can provide full path or place .nii file in supplementary_data

direction = "atlas_to_image";                % "atlas_to_image" or "image_to_atlas". Direction to perform registration. Note for image_to_atlas, inverse must be calculated
calculate_inverse = "false";                 % true, false. Whether to calculate the inverse transform using elastix's Displacement Magnitude Penalty

mask_coordinates = [];                      % Row start, row end, col start, col end. Create a mask around edges during registration
points_file = [];                           % Name of points file to guide registration
save_registered_image = "true";             % Save a copy of registration results

hemisphere = "left";                        % left, right, whole (which brain hemisphere or whole brain)
orientation = "lateral";                    % lateral, dorsal, ventral. Which orientation is at z=0 (i.e. for hemisphere laying flat on midline, this would be lat

%% Nuclei Detection
count_method = "3dunet";              % 3dunet, hessian. 
measure_local_background = "load";    % true, false; measure local image background

use_mask = "true";                    % true, custom, false; Use mask for cell counting
mask_file = [];                       % Name of custom defined mask file. Otherwise, leave empty
mask_resolution = [];                 % Image resolution of mask file

% 3D-Unet
model_file = '128_model.h5';        % Model file name located in /analysis/3dunet/nuclei/models
chunk_size = [112,112,32];          % Chunk size of unet model
chunk_overlap = [16,16,8];          % Overlap between chunks
trained_resolution = [1.21,1.21,4]; % Resolution at which the model was trained. Only required if resampling
resample_chunks = "false";          % Resample images to match model resolution. This process takes significantly longer
gpu = '0';                          % Cuda visible device index

% Hessian blob detector
min_intensity = 200;                    % Minimum intensity for cell nuclei
nuc_diameter_range = [6,20];            % Range of nuclei diameters (in pixels)

%% Cell-Type Classification
classify_method = "";                       % threhsold, gmm, svm; Cell-type classification method

z_normalization = 'false';                  % Apply z normalization to thresholds
log_outliers = 0;                           % Apply log to z scores above or below this threshold
skip_c1 = "true";                           % Skip classifying 1st (reference/nuclei) channel

% Threshold classification
thresholds = [0.005,0.005];                 % Intensity thresholds 
expression{1} = "1*mode + 5*mad";           % Threshold expression
expression{2} = "1*mode + 2*mad";           % Threshold expression

% Gaussian-Mixture Model (gmm) classification
n_clusters = 4;                             % Number of clusters for GMM
mix_proportions = [0.55,0.15,0.30,0.01];	% Intial mixing proportions for n-1 markers for intializing GMM. Leave empty to estimate with k-means++
mix_markers = {[],2,3,[2,3]};               % Marker numbers for the mixing proprtions described above
confidence = [0.5,0.5,0.5,0.5];             % Posterior probability of cell being positive for a marker. This will get adjusted with addition of more clusters
stratify_structures = "true";               % Run GMM clustering on individual structures
split_markers = "false";                    % Run GMM seperately for each marker

% Support-Vector Machine (svm) classification



%% Additional carry-over threshold values from NM_process
lowerThresh = [];                       % 1xn_channels numeric; Lower intensity for rescaling
signalThresh = [];                      % 1xn_channels numeric; Rough estimate for minimal intensity for features of interest
upperThresh = [];                       % 1xn_channels numeric; Upper intensity for rescaling
