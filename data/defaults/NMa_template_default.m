%% Template to run Tissue Clearing Analysis Pipeline
% These are the key parameters
% Set flags to indicate how/whether to run process
resample_images = "true";                       % true, update, false; Perform image resampling 
register_images = "true";                      % true, update, false; Register image to reference atlas.
count_nuclei = "true";                         % true, update, false; Count cell nuclei or other blob objects.
classify_cells = "false";                       % true, update, false; Classify cell-types for detected nuclei centroids.

use_processed_images = "stitched";              % false or name of sub-directory in output directory (i.e. aligned, stitched...); Direct pipeline to load previously processed images in output directory
save_samples = "true";                          % true, false

%% Annotation Parameters
% Orientation key: anterior(a)/posterior(p), superior(s)/inferior(i), left(l)/right(r)
use_annotation_mask = "true";            % true, false or name of file; Use annotation mask for cell counting. Specify mask filename and place in /data/usr/masks/
structures = "cortex";                   % Specify csv file in /annotations detailing which structures to analyze

% For custom annotations
%annotation_resolution = 25;              % Isotropic image resolution of mask file if using custom
%annotation_orientation = 'lps';          % Orientation of custom annotation
%annotation_hemisphere = 'left';          % Hemisphere of annotation

%% Resampling Parameters
resample_resolution = 25;                   % Isotropic resample resolution. This is also the resolution at which registration is performed
resample_channels = 1;                      % Resample specific channels. Only resampled reference channel (nuclei) required for registration
normalize_orientation = "true";             % true, false; Adjust orientation to default atlas orientation

%% Registration Parameters
registration_parameters = [];                   % Name of folder containing elastix registration parameters. Place in /data/elastix_parameter_files/atlas_registration. If empty, will load from "default" or "points"
register_channels = 1;                          % Which channel to register to atlas. Can select more than 1
atlas_file = "ara_nissl_25.nii";                 % Name of the atlas file to register to. Can provide full path or place .nii file in /data/atlas
direction = "atlas_to_image";                   % "atlas_to_image","image_to_atlas","mri_to_image","image_to_mri","mri_to_atlas","atlas_to_mri". Direction to perform registration. Note for image_to_atlas, inverse must be calculated
calculate_inverse = "false";                    % true, false. Whether to calculate the inverse transform from the given registration direction
save_registered_images = "true";                % Whether to save registered images

mask_coordinates = [];                      % Row start, row end, col start, col end. Create a mask around edges during registration
points_file = [];                           % Name of points file to guide registration

%% Nuclei Detection
count_method = "3dunet";                % 3dunet, hessian. 
min_intensity = 200;                    % Minimum intensity for cell nuclei

% 3D-Unet
model_file = '121_model.h5';        % Model file name located in /analysis/3dunet/nuclei/models
chunk_size = [112,112,32];          % Chunk size of unet model
chunk_overlap = [16,16,8];          % Overlap between chunks
trained_resolution = [1.21,1.21,4]; % Resolution at which the model was trained. Only required if resampling
resample_chunks = "false";          % Resample images to match model resolution. This process takes significantly longer
gpu = '0';                          % Cuda visible device index

% Hessian blob detector
nuc_diameter_range = [6,20];            % Range of nuclei diameters (in pixels)

%% Cell-Type Classification
classify_method = "svm";                    % threhsold, gmm, svm; Cell-type classification method
classify_channels = [1,2,3];                % which channels to use for classification
remeasure_centroids = "false";              % true,false; Re-measure channel intensities and save into centroids sheet

% Threshold classification
thresholds = [0.005,0.005];                 % Intensity thresholds 
expression = "1*mode + 5*mad";              % Threshold expression

% Gaussian-Mixture Model (gmm) classification
n_clusters = [];                            % Number of clusters for GMM
split_markers = "true";                     % Run GMM seperately for each marker
mix_proportions = [0.55,0.15,0.30,0.01];	% Intial mixing proportions for n-1 markers for intializing GMM. Leave empty to estimate with k-means++
mix_markers = {[],2,3,[2,3]};               % Marker numbers for the mixing proprtions described above
confidence = [0.5,0.5,0.5,0.5];             % Posterior probability of cell being positive for a marker. This will get adjusted with addition of more clusters
stratify_gmm_csv = "true";                  % Run GMM clustering sper individual structures
%z_normalization = 'false';                  % Apply z normalization to thresholds
%log_outliers = 0;                           % Apply log to z scores above or below this threshold
%skip_c1 = "true";                           % Skip classifying 1st (reference/nuclei) channel

% Support-Vector Machine (svm) classification
load_patches = "true";                          % Load previous centroid patches
load_groups = [];                                % Merge patch annotation from different groups for training
keep_classes = 1:3;                             % Which classes to keep after prediction. Other labeled classes will be discarded
classify_by_annotations = "false";              % Use structure annotations during model training/predictions. Treated as factors
patch_size = [50,6];                            % 1x2 integer. [Patch size for viewing, patch size read by classifier]
n_patches = 1000;                               % integer. Number of patches to generate
min_class_thresh = 0.5;                         % numeric <1. Remove cells dim in every channel from classifier. (i.e. 0.5 removes cells below 50th percentile for all channels)

%% Additional carry-over threshold values from NM_process
lowerThresh = [];                       % 1xn_channels numeric; Lower intensity for rescaling
signalThresh = [];                      % 1xn_channels numeric; Rough estimate for minimal intensity for features of interest
upperThresh = [];                       % 1xn_channels numeric; Upper intensity for rescaling
