%% Template to run NuMorph analysis steps
% Set flags as logical to indicate whether to run process
use_processed_images = "stitched";         % Name of processed image directory (e.g. aligned, stitched). Otherwise set to false

% Image resampling
resample_image = "true";        % true, false, load

% Registration
register_image = "false";        % true, false, load
generate_mask = "false";         % true, false, load

% Nuclei detection
count_cells = "load";                 % 3dunet, hessian, load
use_mask = "true";                    % true, false; Use mask for cell counting
measure_local_background = "load";    % true, false; measure local image background

% Cell-type classification
count_colocalized = "false";          % gmm, svm, threshold, false; Cell-type classification method

save_counts = "true";             % true, false, overwrite 

%% Resampling Parameters
resample_markers = 1:3;        % Specify which channels to resample. Only reference channel 1 needed for registration
resample_res = [10,10,10];     % Resolution to resample to. Should match resolution of ARA reference

%% Registration Parameters
registration_method = "p";      % Affine, BSpline, Points, Other
structures_of_interest = "cortex";       % Structure of interest
atlas_file = 'ara_nissl_left_10.nii';    % Name of the atlas file to register to
hemisphere = "left";            % left, right, whole (which brain hemisphere or whole brain)
orientation = 'lateral';        % lateral, dorsal, ventral. Which position is at z=0 
direction = 'img_to_atlas';         % Forward = register image to atlas; Reverse = register atlas to image
mask_coordinates = [];          % Row start, row end, col start, col end
save_registered_image = 'true';     % Save a copy of registration results

%% Centroid Prediction - 3DUnet
model_file = '128_model.h5';    % Model file name located in /analysis/3dunet/nuclei/models
chunk_size = [112,112,32];      % Chunk size of unet model
chunk_overlap = [16,16,8];           % Overlap between chunks

trained_resolution = [1.21,1.21,4]; % Resolution at which the model was trained. Only required if resampling
resample_chunks = "false";          % Resample images to match model resolution. This process takes significantly longer

gpu = '0';                      % Cuda visible device index

%% Centroid Prediction - Hessian Method
min_intensity = 200;                    % Minimum intensity for cell nuclei
nuc_diameter_range = [6,20];            % Range of nuclei diameters (in pixels)

%% Local Background Measurement
back_res = [25,25,25];

%% Cell-Type Classification
z_normalization = 'false';      % Apply z normalization to thresholds
log_outliers = 0;               % Apply log to z scores above or below this threshold
skip_c1 = "true";               % Skip classifying 1st (reference/nuclei) channel
subtract_background = "false";  % Subtract local background

% Threshold classification
thresholds = [0.005,0.005];           % Intensity thresholds 
expression{1} = "1*mode + 5*mad";    % Threshold expression
expression{2} = "1*mode + 2*mad";    % Threshold expression

% GMM classification
n_clusters = 4;				% Number of clusters for GMM
mix_proportions = [0.55,0.15,0.30,0.01];	% Intial mixing proportions for n-1 markers for intializing GMM. Leave empty to estimate with k-means++
mix_markers = {[],2,3,[2,3]};		% Marker numbers for the mixing proprtions described above
confidence = [0.5,0.5,0.5,0.5];		% Posterior probability of cell being positive for a marker. This will get adjusted with addition of more clusters
stratify_structures = "true";       % Run GMM clustering on individual structures
split_markers = "false";             % Run GMM seperately for each marker

% SVM classification




%% Do not edit these remaining lines
%--------------------------------------------------------------------------
% Save variables and run
temp_path = fileparts(which('NMa_template'));
addpath(genpath(fullfile(temp_path,'..')))
home_path = fileparts(which('NM_analyze.m'));
cd(home_path)
save(fullfile('templates', 'NM_variables.mat'),'-mat')

% Load and append sample info
if exist('sample','var') == 1
    [img_directory, output_directory] = NMsamples(sample);
else
    clear
    error("Sample information is unspecified. Set 'sample' variable.")
end

% Update image directory if using processed images
if ~isequal(use_processed_images,"false")
    img_directory = fullfile(output_directory,use_processed_images);
    if ~exist(img_directory,'dir')
        error("Could not locate processed image directory %s\n",img_directory)
    else
        save(fullfile('templates','NM_variables.mat'),'img_directory','-mat','-append')
    end
end

% Make an output directory
if exist(output_directory,'dir') ~= 7
    mkdir(output_directory);
end

% Make a variables directory
if exist(fullfile(output_directory,'variables'),'dir') ~= 7
    mkdir(fullfile(output_directory,'variables'))
end

% Run analysis
clear
if exist('run_analysis','var') == 1 && run_analysis
    NM_analyze
else
    load 'NM_variables.mat'
    config = load(fullfile('templates','NM_variables.mat'));
end
