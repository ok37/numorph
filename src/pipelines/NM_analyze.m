function [config, path_table] = NM_analyze(config, step)
%--------------------------------------------------------------------------
% NuMorph analysis pipeline designed to perform image registration, nuclei
% counting, and cell-type classification based on nuclear protein markers
% in whole brain images.
%--------------------------------------------------------------------------
% Usage:  
% [config, path_table] = NM_analyze(config, step)
% 
% Inputs:
%   config - config structure for processing
%   step - 'stitch', 'align', 'intensity'. Perform only one of the 3 steps
%   use_adjustments - apply intensity adjustments if performing 1 step
%
% Output:
%   config - config structure after processing
%   path_table - table of filenames used and additional image information
%--------------------------------------------------------------------------

% Load configuration from .mat file, if not provided
if nargin<1
    config = load(fullfile('templates','NM_variables.mat'));
elseif isstring(config)
    config = load(fullfile(config,'NM_variables.mat'));
elseif ~isstruct(config)
    error("Invalid configuartion input")
end
config.var_directory = fullfile(config.output_directory,'variables');

% Default to run full pipeline
if nargin<2
    step = 'analyze';
end

% Make an output directory
if exist(config.output_directory,'dir') ~= 7
    mkdir(config.output_directory);
end

% Make a variables directory
if exist(config.var_directory,'dir') ~= 7
    mkdir(config.var_directory);
end

fprintf("%s\t Working on sample %s \n",datetime('now'),config.sample_name)

%% Create directories
% Update image directory if using processed images
if ~isequal(config.use_processed_images,"false")
    config.img_directory = fullfile(config.output_directory,config.use_processed_images);
    if ~exist(config.img_directory,'dir')
        error("Could not locate processed image directory %s\n",config.img_directory)
    end
end

%% Read image filename information
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate table containing image information
path_table = path_to_table(config);

% Count number of x,y tiles for each channel
nchannels = length(unique(path_table.channel_num));
ntiles = length(unique([path_table.x,path_table.y]));
assert(ntiles == 1, "To perform analysis, there should be only 1 tile for "+...
    "each channel in the image dataset.")

% Check if all resolutions are equal
equal_res = all(cellfun(@(s) config.resolution{1}(3) == s(1,3),config.resolution));

%% Run single step and return if specified
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Note no specific checks on tiles/markers included
if nargin>1 && isequal(step,'resample') 
    config = perform_resampling(config, path_table);
    if nargout == 1; path_table = []; end
    if nargout < 1; clear config; end
    return
elseif nargin>1 && isequal(step,'register')
    config = perform_registration(config);
    if nargout == 1; path_table = []; end
    if nargout < 1; clear config; end
    return
elseif nargin>1 && isequal(step,'count')
    [config, path_table] = perform_channel_alignment(config, path_table, equal_res);
    if nargout == 1; path_table = []; end
    if nargout < 1; clear config; end
    return
elseif nargin>1 && isequal(step,'classify')
    [config, path_table] = perform_channel_alignment(config, path_table, equal_res);
    if nargout == 1; path_table = []; end
    if nargout < 1; clear config; end
    return
end

%% Function for Generating Mask
switch generate_mask
    case 'true'
        fprintf('%s\t Generating mask for selected structures \n',datetime('now'))
        I_mask = gen_mask(hemisphere, structures_of_interest);
        
        if ~isempty(reg_params) && ~isempty(reg_params.atlas_to_img)
            fprintf('%s\t Applying registration parameters \n',datetime('now'))
            
        % Adjust sizes and spacing
        s = 0.4;    % Default: brain registration performed at 25um, mask at 10um
        size1 = reg_params.native_img_size;
        for j = 1:length(reg_params.atlas_to_img.TransformParameters)
            reg_params.atlas_to_img.TransformParameters{j}.FinalBSplineInterpolationOrder = 0;
            reg_params.atlas_to_img.TransformParameters{j}.Size = size1;
            reg_params.atlas_to_img.TransformParameters{j}.Spacing = [s, s, s];
        end
        
        % Transform atlas
        I_mask = transformix(I_mask,reg_params.atlas_to_img,[s, s, s], []);
        
        %%%%%% Additional option for comparing with true positive trace and
        %%%%%% measuring DICE score
        measure_dice = false;
        if measure_dice
            % Find true positive, manually traced mask. This is at 10um resolution
            mask_location = fullfile(output_directory,'resampled');
            files = dir(mask_location);
            file_idx = arrayfun(@(s) contains(s.name,'_mask.tif'), files);

            % Load manually traced mask and resize
            I_true = loadtiff(fullfile(files(file_idx).folder,files(file_idx).name));
            I_true = imresize3(I_true,0.4/s,'Method','nearest');
            
            % Display DICE
            dices = dice(logical(I_true),logical(I_mask));
        end
        
        % Save a .tif copy in the resampled directory
        I_mask_small = imresize3(uint16(I_mask),s,'Method','nearest');
        options.overwrite = true;
        saveastiff(I_mask_small,char(fullfile(output_directory,'resampled','I_mask.tif')),...
            options)
        
        % Save mask as variable
        save(fullfile(output_directory,'variables','I_mask.mat'),'I_mask')
            
        else
            warning('%s\t No registration parameters loaded. Saving mask '+...
                'without applying transformation\n',string(datetime('now')))
        end
    case 'load'
        % Attempt to load annotation mask
        fprintf('%s\t Loading registration parameters \n',datetime('now'))
        try
            load(fullfile(output_directory,'variables','I_mask.mat'),'I_mask')
        catch ME
            error('%s\t Could not locate annotation mask \n',string(datetime('now')))
        end
    case 'false'
        % Use whole image?
        fprintf("%s\t No image mask selected. \n",datetime('now'));
        I_mask = [];
end

%% Read File Information from Stitched Images
% Load images from stitched directory
path_cell = {};
if exist(img_directory,'dir') == 7
    fprintf('%s\t Reading image filename information from stitched directory \n',datetime('now'))
    path_cell{1} = dir(img_directory);
    if isequal(fullfile(output_directory,'stitched'),img_directory)
        location = 'stitched';
        path_table = path_to_table(path_cell,location,[],[]);
    end
else
    error('%s\t Could not locate images in the specified image directory',string(datetime('now')));
end
    
%% Segment and/or Count Cells
switch count_cells
    case "hessian"
        % Detect nuclei centroids using hessian-based blob detector
        % Depcrecated: needs reconfiguring
        path_sub = path_table(path_table.markers == markers(1),:); 
        [pixel_idx_list,box_range] = detect_blobs(path_sub,I_mask, nuc_diameter_range, resolution, adj_param);
        save('pixel_idx_list.mat','pixel_idx_list')
    case "3dunet"
        % Detect nuclei centroids using 3dunet (requires python + conda environment)
        predict_centroids_3dunet(config)
    case "load"
        % Load previos centroid list
        % Default location is: output_directory/(sample_name)_centroids.csv
        path_centroids = fullfile(output_directory, sprintf('%s_centroids.csv',sample_name));
        if exist(path_centroids,'file') == 2
            fprintf('%s\t Loading previous centroid list \n',datetime('now'))
            centroids = readmatrix(path_centroids);
        else
            error("Centroid file %s does not exist",path_centroids)
        end
end

%% Classify Cell-Types
% Load previous cell counting results. Should be located in the output directory.
path_centroids = fullfile(output_directory, sprintf('%s_centroids.csv',sample_name));
if exist(path_centroids,'file') == 2
    fprintf('%s\t Loading previous centroid list \n',datetime('now'))
    centroids = readmatrix(path_centroids);
else
    fprintf("%s\t Could not locate %s in the output directory \n",datetime('now'),path_centroids);
    count_colocalized = 'none';
end

if isequal(config.remeasure_centroids,'true')
    centroids = remeasure_centroids(centroids,path_table,config);
    writematrix(centroids,path_centroids);
end

switch count_colocalized
    case 'gmm'        
        [ct, p, gm] = classify_cells_gmm(centroids, config);
        ceentroids(:,8) = ct;

%        path_centroids = fullfile(output_directory, sprintf('%s_centroids.csv',sample_name));
%        writematrix(centroids,path_centroids);
        
%        I_final = create_image_slices(centroids, path_table_stitched, config);
%        I_final = create_image_slice(centroids, path_table_stitched, config);
        
        %if isequal(config.save_counts,'true') || isequal(config.save_counts,'overwrite')
        %    fprintf('%s\t Saving centroid list \n',datetime('now'))
        %    path_centroids = fullfile(output_directory, sprintf('%s_centroids2_new.csv',sample_name));
        %    writematrix(centroids_new,path_centroids);
        %end
    case 'supervised'
        %get_centroid_patches(centroids,path_table_stitched, config,[50,7])
        ct = classify_cells_svm(centroids, path_table, config);
        
    case 'threshold'
        ct = classify_cells_threshold(centroids, config, B, S);
        centroids(:,end+1) = ct;
    case 'none'
        fprintf('%s\t No centroid list so skipping cell classification \n',datetime('now'))
end

fprintf('%s\t Analysis steps completed! \n',datetime('now'))

end


function config = perform_resampling(config, path_table)
% Image resampling

res = config.resample_resolution;
flags = true(1,length(config.markers));
res_path = fullfile(config.output_directory,'resampled');
config.resampled_paths = cell(1,length(config.markers));

% Create resampled directory
if exist(res_path,'dir') ~= 7
    mkdir(res_path)
else
    % Check for images alread resampled
    for i = config.resample_channels
        filename = fullfile(res_path,sprintf('%s_C%d_%s_%d_%d_%d.nii',...
            config.sample_name,i,config.markers(i),res(1),res(2),res(3)));
        if exist(filename,'file') ~= 2
            flags(i) = false;
        else
            config.resampled_paths{i} = filename;
        end
    end

    % If all channels resampled at correct resolution, return
    if ~isequal(config.resample_images,"update")
        if all(flags)
            return
        else
            config.resample_channels = find(~flags);
        end
    end
end

% Perform resampling
for k = config.resample_channels
    fprintf('%s\t Resampling channel %s\n',datetime('now'),config.markers(k))            
    path_sub = path_table(path_table.markers == config.markers(k),:);
    resample_path_table(path_sub, config);
end

end


function config = perform_registration(config)
% Image registration

fprintf('%s\t Performing image registration \n',datetime('now'))

% Get which structures
[~,structures] = fileparts(config.structures);

% 
reg_dir = fullfile(config.output_directory,'registered');
reg_file = fullfile(reg_dir,'reg_params.mat');
mask_dir = fullfile(config.output_directory,'masks');
mask_file = fullfile(mask_dir,sprintf('%s_C1_%s_%d_%s_mask.nii',...
            config.sample_name,config.markers(1),....
            config.resample_resolution(3),structures));
reg_params = [];

% Create mask directory
if exist(reg_dir,'dir') ~= 7
    mkdir(reg_dir)
end

% Create registerd directory
if exist(reg_dir,'dir') ~= 7
    mkdir(reg_dir)
end

% Attempt to load registration parameters
if exist(reg_file,'file') == 7    
    % Load previous parameters if not updating
    if ~isequal(config.register_image,"update")
        fprintf('%s\t Loading previosuly calculated registration parameters \n',datetime('now'))
        load(reg_file,'reg_params')
    end
    
    % Check if mask exists and return if not updating
    if isequal(config.use_mask,"true") && exist(mask_file,'file') == 7 && ~isempty(reg_params)
        fprintf('%s\t Loading previosuly calculated annotation mask \n',datetime('now'))
        config.mask_path = {mask_file};
        return        
    end
end
        
% Get resampled paths
resample_table = path_to_table(config,'resampled');

% Subset channels to register
idx = ismember(resample_table.markers,config.markers(config.register_channels));
resample_table = resample_table(idx,:);
assert(length(unique(resample_table.x)) == 1 && length(unique(resample_table.y)) == 1 &&...
    length(unique(resample_table.z)) == 1, "All resampled images for registration must be "+...
    "at the same resolution")
config.resample_resolution = [resample_table.x(1), resample_table.y(1), resample_table.z(1)];

% Calculate registration parameters
reg_params = register_to_atlas(config, resample_table.file);

% Apply transformation to other channels in the sample
if isequal(config.save_registered_image,'true')
    for i = 2:height(path_table_resampled)
        fprintf('%s\t Applying transform to %s images\n',datetime('now'),markers(i)); 
        resample_table = path_table_resampled.file{i};
        apply_transform_to_resampled(resample_table,reg_params)
    end
end

% Generate mask
    
        
end