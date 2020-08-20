NMa_variables
NMsamples(sample)
NM_analyze

function NMa_savevars(sample)
%% Template to run Tissue Clearing Processing Pipeline
% Set flags as logical to indicate whether to run process
% To only count cells, set these next 3 to false
resample_image = "false";        % true, false, load
register_image = "false";        % true, false, load. Note: load will only load registration parameters. Set generate_mask to true/load to create/load an annotation mask
generate_mask = "false";         % true, false, load
measure_local_background = "load";  % true, false, load
count_cells = "load";              % 3dunet, hessian, load
count_colocalized = "false";          % gmm, threshold

save_counts = "true";             % true, false, overwrite 

%% Resampling Parameters
resample_markers = 1;          % Specify which channels to resample
resample_res = [25,25,25];     % Resolution to resample to

%% Registration Parameters
registration_method = "p";      % Affine, BSpline, Points, Other
structures_of_interest = "cortex";       % Structure of interest
atlas_file = 'ara_nissl_left_10.nii';    % Name of the atlas file to register to
hemisphere = "left";            % left, right, whole (which brain hemisphere or whole brain)
orientation = 'lateral';        % lateral, dorsal, ventral. Which position is at z=0 
direction = 'img_to_atlas';         % Forward = register image to atlas; Reverse = register atlas to image
mask_coordinates = [];          % Row start, row end, col start, col end
save_registered_image = 'true';     % Save a copy of registration results

%% Nuclei Counting Parameters
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
%expression{1} = "1*mode + 5*mad";  % Threshold expression
%expression{2} = "1*mode + 2*mad";  % Threshold expression

% GMM clustering information
n_clusters = 4;				% Number of clusters for GMM
mix_proportions = [0.55,0.15,0.30,0.01];	% Intial mixing proportions for n-1 markers for intializing GMM. Leave empty to estimate with k-means++
mix_markers = {[],2,3,[2,3]};		% Marker numbers for the mixing proprtions described above
confidence = [0.5,0.5,0.5,0.5];		% Posterior probability of cell being positive for a marker. This will get adjusted with addition of more clusters
stratify_structures = "true";       % Run GMM clustering on each structure individually
split_markers = "false";             % Run GMM seperately for each marker

%% Save variables and run
addpath(genpath('..'))
home_path = fileparts(which('NM_analyze.m'));
cd(home_path)
save -mat NMa_variables.mat
end
