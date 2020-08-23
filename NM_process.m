%This is a image pre-processing program designed to perform channel
%alignment, intensity adjustment, stitching, and resampling on
%light-sheet images, specifically those acquired by LaVision
%Ultramicroscope II. Performing these steps will generate all files
%necessary to run image registration and cell counting in the
%following steps.
clear
% Load configuration from .mat file
load(fullfile('templates','NM_variables.mat'));
config = load(fullfile('templates','NM_variables.mat'));

img_directory = config.img_directory;
output_directory = config.output_directory;

fprintf("%s\t Working on sample %s \n",datetime('now'),config.sample_name)

%% Create directories
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make an output directory
if exist(output_directory,'dir') ~= 7
    mkdir(char(output_directory));
end

% Make a variables directory
if exist(fullfile(output_directory,'variables'),'dir') ~= 7
    mkdir(fullfile(output_directory,'variables'))
end

% Update image directory if using processed images
if ~isequal(config.use_processed_images,"false")
    img_directory = fullfile(output_directory,config.use_processed_images);
    config.img_directory = img_directory;
    if ~exist(img_directory,'dir')
        error("Counld not locate processed image directory %s\n",img_directory)
    end
    if isequal(config.use_processed_images,"aligned")
        fprintf("%s\t Images already aligned. Skipping channel alignment \n",datetime('now'))
        config.channel_alignment = "false";
        channel_alignment = "false";
    end
end

%% Read image filename information
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isequal(config.use_processed_images,"true")
    % Default read images from aligned directory
    img_directory = fullfile(output_directory,"aligned");
    config.channel_alignment = "false";
    channel_alignment = "false";
end

path_cell = {};
if isequal(img_directory,fullfile(output_directory,'aligned'))
    % Start from after multi-channel alignment    
    fprintf("%s\t Reading image filename information from aligned directory \n",datetime('now'))
    path_cell{1} = dir(img_directory);
    location = "aligned";
    path_table = path_to_table(path_cell,location,markers,channel_num,sample_name);
elseif isequal(img_directory, fullfile(output_directory,'stitched'))
    % Start from after stitching
    fprintf("%s\t Reading image filename information from stitched directory \n",datetime('now'))
    path_cell{1} = dir(img_directory);
    location = "stitched";
    path_table = path_to_table(path_cell,location,markers,channel_num,sample_name);
else
    %Start from raw image directory
    fprintf("%s\t Reading image filename information from raw image directory \n",datetime('now'))
    path_cell = cell(1,length(img_directory));
    for i = 1:length(img_directory)
       path_cell{i} = dir(char(img_directory(i)));
    end
    location = "raw";
    path_table = path_to_table(path_cell,location,markers,channel_num,sample_name);
end

% Count number of x,y tiles
ncols = length(unique(path_table.x));
nrows = length(unique(path_table.y));
x_start = min(path_table.x);
y_start = min(path_table.y);
nb_tiles = ncols * nrows;
position_mat = reshape(1:nb_tiles,[nrows,ncols])';

%% Measure Intensity of All Images
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isequal(config.adjust_intensity,"true") || isequal(config.adjust_intensity,'update')
    if isequal(config.adjust_intensity,"true")
        % Define intensity adjustment parameters from scratch. Intensity
        % adjustment measurements should be made on raw images.
        fprintf("%s\t Defining new adjustment parameters \n",datetime('now'));        
        adj_fields = {'adjust_intensity','adjust_ls_width','adjust_tile_intensity',...
        'rescale_intensities','shading_correction','lowerThresh','upperThresh','gamma'};

        % Update parameters from config structure
        for i = 1:length(adj_fields)
            try
                adj_params.(adj_fields{i}) = config.(adj_fields{i});
            catch
            end
        end
    else
        % Load previous adjustment parameter structure
        fprintf("%s\t Updating previous adjustment parameter fields \n",datetime('now'));
        load(fullfile(output_directory,'variables','adj_params.mat'))
        
        % Update which adjustment to apply based on current configs
        adj_params.adjust_ls_width = config.adjust_ls_width;
        adj_params.adjust_tile_intensity = config.adjust_tile_intensity;
        adj_params.rescale_intensities = config.rescale_intensities;
        adj_params.shading_correction = config.shading_correction;
    end

    % Measure sample images in overlapping regions and calculate
    % intensity thresholds
    lowerThresh_measured = zeros(1,3); upperThresh_measured = zeros(1,3);
    flatfield = cell(1,3);  darkfield = cell(1,3);
    for k = 1:length(markers)            
        fprintf("%s\t Measuring Intensity for %s \n",datetime('now'),config.markers(k));

        stack = path_table(path_table.markers == config.markers(k),:);

        [t_adj, lowerThresh_measured(k), upperThresh_measured(k), y_adj, flatfield{k}, darkfield{k}] = ...
            measure_images(stack, config, k);

        % Update adjustments parameters structure
        if isequal(config.adjust_ls_width,"true") || isequal(config.adjust_intensity,"true")
            adj_params.y_adj{k} = y_adj;
        end
        if isequal(config.adjust_tile_intensity,"true") || isequal(config.adjust_intensity,"true")
            adj_params.t_adj{k} = t_adj;
        end
        if isequal(config.shading_correction,"true") || isequal(config.adjust_intensity,"true")
            adj_params.flatfield{k} = flatfield{k};
            adj_params.darkfield{k} = darkfield{k};
        end
    end

    % Check for user-defined intensity threshold values
    [adj_params, config, lowerThresh, upperThresh] = check_adj_parameters(adj_params,...
        config, lowerThresh_measured, upperThresh_measured);
    config.adj_params = adj_params;

    % Save adjustment parameters to output directory
    fprintf("%s\t Saving adjustment parameters \n",datetime('now'));
    save(fullfile(output_directory,'variables','adj_params.mat'), 'adj_params')
        
    % Save flatfield and darkfield as seperate variables
    if isequal(config.shading_correction,"true")
        fprintf("%s\t Saving flatfield and darkfield images \n",datetime('now'));
        save(fullfile(output_directory,'variables','flatfield.mat'), 'flatfield')
        save(fullfile(output_directory,'variables','darkfield.mat'), 'darkfield')
    end

    config.adjust_intensity = "true";

elseif isequal(config.adjust_intensity,"load")
    % Load adjustment parameters from output directory and use as is  
    fprintf("%s\t Loading adjustment parameters \n",datetime('now'));
    load(fullfile(output_directory,'variables','adj_params.mat'))
    [adj_params, config] = check_adj_parameters(adj_params,config);
    config.adj_params = adj_params;
    config.adjust_intensity = "true";
        
    % Update which adjustment to apply based on current configs
    config.adj_params.adjust_ls_width = config.adjust_ls_width;
    config.adj_params.adjust_tile_intensity = config.adjust_tile_intensity;
    config.adj_params.rescale_intensities = config.rescale_intensities;
    config.adj_params.shading_correction = config.shading_correction;
        
elseif isequal(config.adjust_intensity,"false")
    fprintf("%s\t No intensity Adjustments Selected \n",datetime('now'));
    config.adj_params = [];
        
    if ~isempty(config.lowerThresh) && ~isempty(config.upperThresh)
        config.lowerThresh = config.lowerThresh/65535;
        config.upperThresh = config.upperThresh/65535;
    else
        fprintf("%s\t Values for lower and upper thresholds must be defined. "+...
            "Attempting to load from adjustment parameters.\n",datetime('now'));
        try
            load(fullfile(output_directory,'variables','adj_params.mat'))
        catch ME
            
        end
    end
else
    error("%s\t Unrecognized selection for adjust_intensity. "+...
        "Please select ""true"", ""load"", or ""false"".\n",string(datetime('now')))
end

%% Channel Alignment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch channel_alignment
    case 'translation'
        fprintf("%s\t Aligning channels by translation \n",datetime('now'));
        % Determine z displacement from reference channel
        % Check first if .mat file exists in output directory
        if exist(fullfile(output_directory,'z_displacement_align.mat'),'file') == 2
            fprintf("%s\t Loading z displacement matrix \n",datetime('now'));
            load(fullfile(output_directory,'z_displacement_align.mat'))
        else
            z_displacement_align = zeros([nrows,ncols,numel(markers)-1]);
            for k = 2:length(markers)
                fprintf("%s\t Performing channel alignment z calculation for %s/%s \n",...
                    datetime('now'),markers{k},markers{1});
                z_tile = zeros(1,numel(position_mat));
                ave_signal = zeros(1,numel(position_mat));
                for idx = 1:numel(position_mat)
                    [y,x] = find(position_mat==idx);
                    path_ref = path_table(path_table.x==x & path_table.y==y & path_table.markers == markers{1},:);
                    path_mov = path_table(path_table.x==x & path_table.y==y & path_table.markers == markers{k},:);
                    % This measures the displacement in Z for a given channel to the reference
                    [z_tile(idx),ave_signal(idx)] = zAlign2(path_mov,path_ref,z_positions,z_window,...
                        lowerThresh(k),z_initial(k));
                    
                    fprintf("%s\t Predicted z displacement of %d for position %d \n",...
                        datetime('now'),z_tile(idx),idx);
                end
                z_matrix = reshape(z_tile,[nrows, ncols])';
                z_displacement_align(:,:,k-1) = z_matrix;
            end
            % Save displacement variable to output directory
            save(fullfile(output_directory,'z_displacement_align.mat'), 'z_displacement_align')
        end

        % Perform channel alignment
        alignment_table = cell(ncols,nrows);
        for idx = 1:numel(position_mat)
            [y,x] = find(position_mat==idx);            
            fprintf("%s\t Aligning channels to %s for tile 0%dx0%d \n",...
                        datetime('now'),markers{1},y,x);
            path_align = path_table(path_table.x==x & path_table.y==y,:);
            alignment_table{y,x} = coregister(path_align,z_displacement_align(y,x,:),...
                output_directory,lowerThresh,align_stepsize,save_aligned_images);
        end
        
        % Save registration variable to output directory
        save(fullfile(output_directory,'alignment_table.mat'), 'alignment_table')
        fprintf("%s\t Alignment completed! \n",datetime('now'));
        
        % Change image directory to aligned directory so that subsequent
        % steps load these images
        if isequal(save_aligned_images,"true")
            img_directory = fullfile(char(output_directory),'aligned');
        end
    case 'elastix'
        fprintf("%s\t Aligning channels using B-splines \n",datetime('now'));

        % Define which tiles to align if only doing subset
        if isempty(config.align_tiles)
            tiles_to_align = 1:nb_tiles;
        else
            tiles_to_align = config.align_tiles;
        end
        
        % Check for intensity bounds which are required for alignment
        if isequal(config.load_alignment_params,"false") && isempty(config.lowerThresh) || isempty(config.upperThresh)
            fprintf("%s\t Lower and upper intensity thresholds are unspecified but are required "+...
                "for accurate alignment. Measuring these now... \n",datetime('now'));
            for k = 1:length(markers)            
                fprintf("%s\t Measuring Intensity for %s \n",datetime('now'),config.markers(k));
                stack = path_table(path_table.markers == config.markers(k),:);
                adj_params = [];
                [lowerThresh(k), upperThresh(k)] = get_thresholds(stack,config);
                fprintf("%s\t Calculated lower threshold: %.2f \n",datetime('now'),lowerThresh(k));
                fprintf("%s\t Calculated upper threshold: %.2f \n",datetime('now'),upperThresh(k));
            end
            config.lowerThresh = lowerThresh/65535;
            config.upperThresh = upperThresh/65535;            
        end

        % Perform alignment
        for i = tiles_to_align
            [y,x] = find(position_mat==i);
            path_align = path_table(path_table.x==x & path_table.y==y,:);
            fprintf("%s\t Aligning channels to %s for tile 0%dx0%d \n",...
                datetime('now'),markers{1},y,x);
            elastix_channel_alignment(x, y, path_align, config);
        end
                
        % Change image directory to aligned directory so that subsequent
        % steps load these images
        config.img_directory = fullfile(output_directory,'aligned');
        img_directory = fullfile(output_directory,'aligned');
        path_cell = {dir(config.img_directory)};
        location = "aligned";
        path_table = path_to_table(path_cell,location,markers(1),channel_num,sample_name);
        
        % Update tile intensity adjustments using newly aligned images.
        % Also, set light sheet width adjustments + flatfield adjustments
        % to false as these were applied during the alignment step
        if ~isempty(config.adj_params)
            config.adjust_ls_width = "false";
            config.adj_params.adjust_ls_width = "false";
            config.shading_correction = "false";
            config.adj_params.shading_correction = "false";

            for k = 1:length(markers)            
                fprintf("%s\t Updating intensity measurements for marker %s using "+...
                    "newly aligned images \n",datetime('now'),config.markers(k));

                stack = path_table(path_table.markers == config.markers(k),:);

                t_adj = measure_images(stack, config, k);

                adj_params.t_adj{k} = t_adj;
            end
        end
        
    case "false"
        % No channel alignment
        fprintf("%s\t No channel alignment selected \n",datetime('now'));
end

%% Stitching
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isequal(stitch_img,"false")
    fprintf("%s\t No image stitching selected \n",datetime('now'));
else
    fprintf("%s\t No image stitching selected \n",datetime('now'));
    % Check for aligned images
    if isequal(config.img_directory, fullfile(output_directory,'aligned'))
        fprintf("%s\t Using images from aligned directory \n",datetime('now'));
        path_reg{1} = dir(img_directory);
        location = 'aligned';
        path_table = path_to_table(path_reg,location,markers(1),channel_num,sample_name);

        %%%%%%%%Edit this 
    elseif length(config.markers) > 1 && exist(fullfile(config.output_directory,'alignment_table.mat'),'file') == 2 && ~exist('alignment_table','var')
        fprintf("%s\t Using multi-channel alignment measurements from .mat file \n",datetime('now'));
        load(fullfile(char(output_directory),'alignment_table.mat'))
        path_table = cellfun(@(s) s{:,2}, alignment_table,'UniformOutput',false);
        alignment_table = cellfun(@(s) s{:,1}, alignment_table,'UniformOutput',false);
        path_table = vertcat(path_table{:});
    end
  
    %Trim z positions that aren't present in all tiles
    min_z_mat = zeros(nrows,ncols);
    max_z_mat = zeros(nrows,ncols);
    for i = 1:ncols
        for j = 1:nrows
            max_z_mat(j,i) = max(path_table.z(path_table.x==i & path_table.y==j));
            min_z_mat(j,i) = min(path_table.z(path_table.x==i & path_table.y==j));
        end
    end
    max_z = min(max_z_mat(:));
    min_z = max(min_z_mat(:));
    path_table(path_table.z>max_z,:)=[];
    path_table(path_table.z<min_z,:)=[];

    % Adjust z position based overlapping regions 
    switch adjust_z 
        case "true"
            % Calculate z displacements
            fprintf("%s\t Calculating z displacements for stitching \n",datetime('now'));
            z_adj = calculate_adjusted_z(path_table,nrows,ncols,markers,...
                overlap,z_positions,z_window,lowerThresh,output_directory);
            [~,z_idx] = setdiff(path_table.file,z_adj.file);
            path_table(z_idx,:) = [];
            path_table.z_adj = z_adj.z_adj;

        case 'file'
            % Load previously calculated z displacements for these images 
            fprintf("%s\t Loading adjusted z positions \n",datetime('now'));
            load(fullfile(output_directory,'variables','adjusted_z.mat'), 'z_adj')
            [~,z_idx] = setdiff(path_table.file,z_adj.file);
            path_table(z_idx,:) = []; 
            try
                path_table.z_adj = z_adj.z_adj;
            catch ME
                if isequal(ME.identifier,'MATLAB:table:RowDimensionMismatch')
                    error("%s\t Loaded z-positions do not match specified "+...
                        "file information. Rerun adjust_z \n",string(datetime('now')))
                end
            end
        case 'matrix'
            % Load z displacement info from a matrix, ignoring whether
            % filenames match up
            fprintf("%s\t Loading z displacement matrix \n",datetime('now'));
            load(fullfile(output_directory,'variables', 'z_disp_matrix'), 'z_disp_matrix')
            z_adj = apply_adjusted_z(path_table, z_disp_matrix);
            [~,z_idx] = setdiff(path_table.file,z_adj.file);
            path_table(z_idx,:) = []; 
            path_table.z_adj = z_adj.z_adj;
        case "false"
            path_table.z_adj = path_table.z;
    end
    
    % Trim markers to only the ones specified
    if ~isempty(config.stitch_sub_channel)
       markers = markers(config.stitch_sub_channel);
       path_table = path_table(ismember(path_table.markers,markers),:);
    end
    
    % Create empty alignment_table if not doing multi-channel alignment
    if exist('alignment_table','var') == 0
        config.alignment_table = [];
    end
        
    % Check whether to stitch from previously calculated transforms.
    % Otherwise stitch from scratch
    if isequal(stitch_img,"load")
        fprintf("%s\t Loading previously calculated stitching parameters \n",datetime('now'));

        % Load stitching parameters
        load(fullfile(output_directory,'variables','stitch_tforms.mat'))
                
        % Check if sizes match up
        nb_images = length(min(path_table.z_adj):max(path_table.z_adj));

        if size(h_stitch_tforms,2) ~= nb_images
           error("%s\t Number of slices z slices in stitching parameters does not match loaded images \n",string(datetime('now')))
        elseif size(h_stitch_tforms,1) ~= (ncols-1)*nrows*2
           error("%s\t Number of columns positions does not match loaded stitching parameters \n",string(datetime('now')))
        elseif size(v_stitch_tforms,1) ~= (nrows-1)*2
           error("%s\t Number of row positions does not match loaded stitching parameters \n",string(datetime('now')))
        end
        
        % Perform stitching using loaded parameters
        stitch_from_loaded_parameters(path_table, h_stitch_tforms, v_stitch_tforms, config)
    else
        % Perform stitching
        stitch_old4(path_table,output_directory,overlap,sd,config,adj_params,...
            markers,stitch_sub_stack,alignment_table,use_middle,subtract_img_background,...
            sift_refinement,blending_method,number_of_cores)
    end

    fprintf("%s\t Stitching completed! \n",datetime('now'));
end
