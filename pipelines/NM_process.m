function NM_process(config)
%--------------------------------------------------------------------------
% NuMorph processing pipeline designed to perform channel alignment,
% intensity adjustment, and stitching on multi-channel light-sheet images.
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
    if isequal(config.use_processed_images,"aligned")
        warning("%s\t Images already aligned. Skipping channel alignment \n",datetime('now'))
        config.channel_alignment = "false";
    end
end

%% Read image filename information
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isequal(config.img_directory,fullfile(config.output_directory,'aligned'))
    % Start from after multi-channel alignment    
    fprintf("%s\t Reading image filename information from aligned directory \n",datetime('now'))
    location = "aligned";
elseif isequal(config.img_directory, fullfile(config.output_directory,'stitched'))
    % Start from after stitching
    fprintf("%s\t Reading image filename information from stitched directory \n",datetime('now'))
    location = "stitched";
else
    %Start from raw image directory
    fprintf("%s\t Reading image filename information from raw image directory \n",datetime('now'))
    location = "raw";
end
path_table = path_to_table(config,location);

% Count number of x,y tiles
ncols = length(unique(path_table.x));
nrows = length(unique(path_table.y));
nb_tiles = ncols * nrows;
position_mat = reshape(1:nb_tiles,[nrows,ncols])';

%% Measure Intensity of All Images
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch config.adjust_intensity
    case {'true', 'update'}
        if isequal(config.adjust_intensity,"true")
            % Define intensity adjustment parameters from scratch. Intensity
            % adjustment measurements should be made on raw images.
            fprintf("%s\t Defining new adjustment parameters \n",datetime('now'));        
            adj_fields = {'adjust_intensity','adjust_tile_shading',...
                'adjust_tile_position','lowerThresh','upperThresh'};
        
            % Update parameters from config structure
            for i = 1:length(adj_fields)
                adj_params.(adj_fields{i}) = config.(adj_fields{i});
            end
        else
            % Load previous adjustment parameter structure
            fprintf("%s\t Updating previous adjustment parameter fields \n",datetime('now'));
            load(fullfile(config.output_directory,'variables','adj_params.mat'),'adj_params')
        end
        
        % Measure sample images in overlapping regions and calculate
        % intensity thresholds for each marker
        lowerThresh_measured = zeros(1,length(config.markers)); upperThresh_measured = lowerThresh_measured;
        flatfield = cell(1,length(config.markers));  darkfield = flatfield;
        for k = 1:length(config.markers)            
            if isequal(config.adjust_intensity,"update") && ~ismember(k,config.update_intensity_channels)
                continue
            end
            fprintf("%s\t Measuring Intensity for %s \n",datetime('now'),config.markers(k));
            stack = path_table(path_table.markers == config.markers(k),:);

            % Take measurements
            [t_adj, lowerThresh_measured(k), upperThresh_measured(k), y_adj, flatfield{k}, darkfield{k}] = ...
                measure_images(config, stack, k);

            % Save new adjustments parameters structure
            if isequal(config.adjust_intensity,"true")
                adj_params.y_adj{k} = y_adj;
                adj_params.t_adj{k} = t_adj;
                adj_params.flatfield{k} = flatfield{k};
                adj_params.darkfield{k} = darkfield{k};
                adj_params.darkfield_intensity{k} = config.darkfield_intensity(k);
            else
                % Update adjustments parameters structure
                if ~all(y_adj==1)
                   adj_params.y_adj{k} = y_adj;
                end
                if ~all(t_adj(:)==1)
                   adj_params.t_adj{k} = t_adj;
                end
                if ~all(flatfield{k}(:)==1)
                    adj_params.flatfield{k} = flatfield{k};
                    adj_params.darkfield{k} = darkfield{k};
                end
                adj_params.darkfield_intensity{k} = config.darkfield_intensity(k);
            end
        end

        % Check for user-defined intensity threshold values
        [adj_params, config, lowerThresh, upperThresh] = check_adj_parameters(adj_params,...
            config, lowerThresh_measured, upperThresh_measured);
        config.adj_params = adj_params;

        % Save adjustment parameters to output directory
        fprintf("%s\t Saving adjustment parameters \n", datetime('now'));
        save(fullfile(config.output_directory,'variables','adj_params.mat'), 'adj_params')

        % Save flatfield and darkfield as seperate variables
        if isequal(config.adjust_tile_shading,"basic")
            fprintf("%s\t Saving flatfield and darkfield images \n",datetime('now'));
            save(fullfile(config.output_directory,'variables','flatfield.mat'), 'flatfield')
            save(fullfile(config.output_directory,'variables','darkfield.mat'), 'darkfield')
        end
        config.adjust_intensity = "true";

    case 'load'
        % Load adjustment parameters from output directory and use as is  
        fprintf("%s\t Loading adjustment parameters \n",datetime('now'));
        load(fullfile(config.output_directory,'variables','adj_params.mat'),'adj_params')
        [adj_params, config] = check_adj_parameters(adj_params,config);
        config.adj_params = adj_params;
        config.adjust_intensity = "true";

        % Update which adjustment to apply based on current configs
        config.adj_params.adjust_tile_shading = config.adjust_tile_shading;
        config.adj_params.adjust_tile_position = config.adjust_tile_position;
        
    case 'false'
        % No intensity adjustments
        fprintf("%s\t No intensity Adjustments Selected \n",datetime('now'));
        config.adj_params = [];

        if ~isempty(config.lowerThresh) && ~isempty(config.upperThresh)
            config.lowerThresh = config.lowerThresh/65535;
            config.upperThresh = config.upperThresh/65535;
        elseif exist(fullfile(config.output_directory,'variables','adj_params.mat'),'file') == 2
            fprintf("%s\t Values for lower and upper thresholds must be defined. "+...
                "Loading these from adjustment parameters that already exist. \n",datetime('now'));
            load(fullfile(config.output_directory,'variables','adj_params.mat'),'adj_params')
            config.lowerThresh = adj_params.lowerThresh;
            config.upperThresh = adj_params.upperThresh;
        end
        
    otherwise
        error("%s\t Unrecognized selection for adjust_intensity. "+...
            "Please select ""true"", ""update"",""load"", or ""false"".\n",string(datetime('now')))
end

%% Channel Alignment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch config.channel_alignment
    case 'translation'
        fprintf("%s\t Aligning channels by translation \n",datetime('now'));
        
        % Determine z displacement from reference channel
        % Check first if .mat file exists in output directory
        if exist(fullfile(config.output_directory,'variables','z_displacement_align.mat'),'file') == 2 || ...
                isequal(config.update_z_adjustment,"true")
            fprintf("%s\t Loading z displacement matrix \n",datetime('now'));
            load(fullfile(config.output_directory,'variables','z_displacement_align.mat'),'z_displacement_align')
        else
            z_displacement_align = zeros([nrows,ncols,numel(config.markers)-1]);
            for k = 2:length(config.markers)
                fprintf("%s\t Performing channel alignment z calculation for %s/%s \n",...
                    datetime('now'),config.markers{k},config.markers{1});
                z_tile = zeros(1,numel(position_mat));
                ave_signal = zeros(1,numel(position_mat));
                for idx = 1:numel(position_mat)
                    [y,x] = find(position_mat==idx);
                    path_ref = path_table(path_table.x==x & path_table.y==y & path_table.markers == config.markers{1},:);
                    path_mov = path_table(path_table.x==x & path_table.y==y & path_table.markers == config.markers{k},:);
                    % This measures the displacement in z for a given channel to the reference
                    [z_tile(idx),ave_signal(idx)] = z_align_channel(path_mov,path_ref,config.z_positions,config.z_window,...
                        config.lowerThresh(k),config.z_initial(k));
                    fprintf("%s\t Predicted z displacement of %d for tile %d x %d \n",...
                        datetime('now'),z_tile(idx),y,x);
                end
                z_matrix = reshape(z_tile,[nrows, ncols])';
                z_displacement_align(:,:,k-1) = z_matrix;
            end
            % Save displacement variable to output directory
            fprintf("%s\t Saving z displacement matrix \n",datetime('now'));
            save(fullfile(config.output_directory,'variables','z_displacement_align.mat'), 'z_displacement_align')
        end

        % Create variable
        alignment_table = cell(ncols,nrows);
        save_path = fullfile(config.output_directory,'variables','alignment_table.mat');
        if ~exist(save_path,'file')
            save(save_path,'alignment_table')
        end
        
        % Define which tiles to align if only doing subset
        tiles_to_align = 1:nb_tiles;
        if ~isempty(config.align_tiles) && all(ismember(config.align_tiles,tiles_to_align))
            tiles_to_align = config.align_tiles;
        elseif ~isempty(config.align_tiles) && ~all(ismember(config.align_tiles,tiles_to_align))
            error("Selected subset of tiles to align outside range of all tiles")
        end
        
        % Perform channel alignment
        for idx = tiles_to_align
            [y,x] = find(position_mat==idx);            
            fprintf("%s\t Aligning channels to %s for tile %d x %d \n",...
                        datetime('now'),config.markers{1},y,x);
            path_align = path_table(path_table.x==x & path_table.y==y,:);
            alignment_table{y,x} = align_by_translation(config,path_align,z_displacement_align(y,x));
            save(save_path,'alignment_table')
        end
        fprintf("%s\t Alignment completed! \n",datetime('now'));
        
        % Change image directory to aligned directory so that subsequent
        % steps load these images
        if isequal(config.save_aligned_images,"true")
            config.img_directory = fullfile(config.output_directory,"aligned");
            location = "aligned";
            path_table = path_to_table(config,location);
        end
    case 'elastix'
        fprintf("%s\t Aligning channels using B-splines \n",datetime('now'));
        
        % Create variable
        alignment_params = cell(ncols,nrows);
        save_path = fullfile(config.output_directory,'variables','alignment_params.mat');
        if ~exist(save_path,'file')
            save(save_path,'alignment_params','-v7.3')
        else
            load(save_path,'alignment_params')
            assert(size(alignment_params,1) == nrows, "Number of row positions do not "+...
                "match loaded alignment parameters")
            assert(size(alignment_params,2) == ncols, "Number of row positions do not "+...
                "match loaded alignment parameters")
        end

        % Define which tiles to align if only doing subset
        tiles_to_align = 1:nb_tiles;
        if ~isempty(config.align_tiles) && all(ismember(config.align_tiles,tiles_to_align))
            tiles_to_align = config.align_tiles;
        elseif ~isempty(config.align_tiles) && ~all(ismember(config.align_tiles,tiles_to_align))
            error("Selected subset of tiles to align outside range of all tiles")
        end
        
        % Check for intensity bounds which are required for alignment
        if isequal(config.load_alignment_params,"false") && isempty(config.lowerThresh) || isempty(config.upperThresh)
            fprintf("%s\t Lower and upper intensity thresholds are unspecified but are required "+...
                "for accurate alignment. Measuring these now... \n",datetime('now'));
            for k = 1:length(config.markers)            
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
                datetime('now'),config.markers{1},y,x);
            alignment_params(y,x) = elastix_channel_alignment(config,path_align,true);
            save(save_path,'alignment_params','-v7.3')
        end
                
        % Change image directory to aligned directory so that subsequent
        % steps load these images
        config.img_directory = fullfile(config.output_directory,'aligned');
        location = "aligned";
        path_table = path_to_table(config,location);
        
        % Update tile intensity adjustments using newly aligned images.
        % Also, set light sheet width adjustments + flatfield adjustments
        % to false as these were applied during the alignment step
        config.adjust_ls_width = "false";
        config.adj_params.adjust_ls_width = "false";
        config.shading_correction = "false";
        config.adj_params.shading_correction = "false";
            
       % In case applying tile adjustments, re-calculate thresholds
        if isequal(config.adjust_intensity,"true") && isequal(config.adjust_tile_intensity,"true")
            for k = 1:length(config.markers)            
                fprintf("%s\t Updating intensity measurements for marker %s using "+...
                    "newly aligned images \n",datetime('now'),config.markers(k));
                stack = path_table(path_table.markers == config.markers(k),:);
                t_adj = measure_images(config, stack, k);
                adj_params.t_adj{k} = t_adj;
            end
        end
    case "false"
        % No channel alignment
        fprintf("%s\t No channel alignment selected \n",datetime('now'));
    otherwise
        error("%s\t Unrecognized selection for channel_alignment. "+...
            "Please select ""translation"", ""elastix"", or ""false"".\n",string(datetime('now')))
end

%% Stitching
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isequal(config.stitch_img,"false")
    fprintf("%s\t No image stitching selected \n",datetime('now'));
else
    fprintf("%s\t Perfoming iterative 2D image stitching \n",datetime('now'));
    
    % Check for aligned images
    if isequal(config.img_directory, fullfile(config.output_directory,'aligned'))
        fprintf("%s\t Using images from aligned directory \n",datetime('now'));
    elseif isequal(config.load_alignment_params,"true")
        fprintf("%s\t Loading channel alignment translations \n",datetime('now'));
        load(fullfile(config.output_directory,'variables','alignment_table.mat'),'alignment_table')
    end
  
    % Trim z positions that aren't present in all tiles
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
    
    % If applying channel alignment, add field and check tiles
    if exist('alignment_table','var') == 0 || isequal(config.img_directory,fullfile(config.output_directory,'aligned'))
        config.alignment_table = [];
    else
        empty_alignment_tiles = find(cellfun(@(s) isempty(s),alignment_table'));
        if empty_alignment_tiles
            error("No alignment parameters to apply for tiles %s",...
                sprintf("%s",num2str(empty_alignment_tiles)))
        end
        config.alignment_table = alignment_table;
    end

    % Calculate adjustments in z
    if isequal(config.update_z_adjustment,"true") || exist(fullfile(config.var_directory,'z_disp_matrix.mat'),'file') ~= 2
        % Calculate z displacements
        fprintf("%s\t Calculating z displacement matrix for stitching \n",datetime('now'));
        z_adj = calculate_adjusted_z(path_table,nrows,ncols,config.markers,...
                config.overlap,config.z_positions,config.z_window,config.lowerThresh,...
                config.output_directory);
        [~,z_idx] = setdiff(path_table.file,z_adj.file);
        path_table(z_idx,:) = [];
        path_table.z_adj = z_adj.z_adj;
    else
        % Load z displacement info from a matrix, ignoring whether
        % filenames match up
        fprintf("%s\t Loading z displacement matrix for stitching \n",datetime('now'));
        load(fullfile(config.var_directory, 'z_disp_matrix.mat'), 'z_disp_matrix')
        assert(nrows == size(z_disp_matrix,1) && ncols == size(z_disp_matrix,2), "Loaded z adjustment parameters do not match number of tiles "+...
            "detected in the input image directory. Check input image directory or update z adjustment")
        z_adj = apply_adjusted_z(path_table, z_disp_matrix);
        [~,z_idx] = setdiff(path_table.file,z_adj.file);
        path_table(z_idx,:) = []; 
        path_table.z_adj = z_adj.z_adj;
    end
    
    % Trim markers to only the ones specified
    if ~isempty(config.stitch_sub_channel)
       markers = config.markers(config.stitch_sub_channel);
       path_table = path_table(ismember(path_table.markers,markers),:);
    end
        
    % Check whether to stitch from previously calculated transforms.
    % Otherwise stitch from scratch
    if isequal(config.stitch_img,"load")
        fprintf("%s\t Loading previously calculated stitching parameters \n",datetime('now'));
        % Load stitching parameters
        try
            load(fullfile(config.var_directory,'stitch_tforms.mat'), 'h_stitch_tforms','v_stitch_tforms')
        catch 
            error("Could not locate stitching transforms in variables directory. Re-run image stitching")
        end
                
        % Check if sizes match up
        nb_images = length(min(path_table.z_adj):max(path_table.z_adj));

        if size(h_stitch_tforms,2) ~= nb_images
           error("%s\t Number of slices z slices in stitching parameters does not match loaded images \n",string(datetime('now')))
        elseif size(h_stitch_tforms,1) ~= (ncols-1)*nrows*2
           error("%s\t Number of column positions does not match loaded stitching parameters \n",string(datetime('now')))
        elseif size(v_stitch_tforms,1) ~= (nrows-1)*2
           error("%s\t Number of row positions does not match loaded stitching parameters \n",string(datetime('now')))
        end
        
        % Perform stitching using loaded parameters
        stitch_from_loaded_parameters(path_table, h_stitch_tforms, v_stitch_tforms, config)
    else
        % Perform stitching
        stitch_iterative(config,path_table)
    end

    fprintf("%s\t Stitching completed! \n",datetime('now'));
end

end