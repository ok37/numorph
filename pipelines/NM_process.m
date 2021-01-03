function [config, path_table] = NM_process(config, step, use_adjustments)
%--------------------------------------------------------------------------
% NM_process NuMorph processing pipeline designed to perform channel 
% alignment, intensity adjustment, and stitching on multi-channel 
% light-sheet images.
%
% Syntax:  [config, path_table] = NM_process(config, step, use_adjustments)
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

% Applying (not calculating) intensity adjustments
if nargin<3
    use_adjustments = true;
elseif use_adjustments
    %config.adjust_intensity = strrep(config.adjust_intensity,"true", "load");
end

% Make an output directory
if exist(config.output_directory,'dir') ~= 7
    mkdir(config.output_directory);
end

% Make a variables directory
config.var_directory = fullfile(config.output_directory,'variables');
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
        warning("Images already aligned. Skipping channel alignment")
        config.channel_alignment = "false";
    end
end

%% Read image filename information
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate table containing image information
path_table = path_to_table(config);

% Count number of x,y tiles for each channel
nchannels = length(unique(path_table.channel_num));
ncols = zeros(1,nchannels); nrows = zeros(1,nchannels);
for i = 1:nchannels
    path_sub = path_table(path_table.channel_num == i,:);
    ncols(i) = length(unique(path_sub.x));
    nrows(i) = length(unique(path_sub.y));
end

% Check if all resolutions are equal
equal_res = all(cellfun(@(s) config.resolution{1}(3) == s(1,3),config.resolution));

%% Intensity Adjustment
% Configure intensity adjustments
if use_adjustments || isequal(step,'intensity')
    [config, path_table] = intensity_adjustment(config, path_table, nrows, ncols);
end

%% Run single step and return if specified
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Note no specific checks on tiles/markers included
if nargin>1 && isequal(step,'stitch')    
    config.stitch_img = 'true';
    [config, path_table] = perform_stitching(config, path_table);
    return
elseif nargin>1 && isequal(step,'align')
    [config, path_table] = align_channels(config, path_table, equal_res);
    return
elseif nargin>1 && isequal(step,'intensity')
    return
end

%% Run multiple steps
% If equal numbers of column and row tiles, perform alignment on each tile,
% then do stitching on aligned images. This assumes the same resolution for
% each channel. If the channels are imaged at different resolutions, stitch
% each channel if multi-tile. Then run alignment
if all(ncols == max(ncols)) && all(nrows == max(nrows))
    % Equal number of tiles
    % Channel Alignment
    [config, path_table] = align_channels(config,path_table,equal_res);
    
    % Stitching
    if ncols(1)*nrows(1) == 1
        fprintf("%s\t 1 tile detected for marker %s \n",datetime('now'),config.markers(1));
    else
        [config, path_table] = perform_stitching(config, path_table);
    end
else
    % Different number of tiles
    fprintf("%s\t Different number of tiles between channels \n",datetime('now'))

    % Check resolutions are present
    assert(length(config.resolution) == length(config.markers),"Specify resolution "+...
        "for each channel as cell array");
    
    % Stitching
    for i = 1:length(config.markers)
        % Continue if only 1 tile
        if ncols(i)*nrows(i) == 1
            fprintf("%s\t 1 tile detected for marker %s \n",datetime('now'),config.markers(i));
            continue
        else
            config.stitch_sub_channel = i;
        end
        [config, path_table] = perform_stitching(config, path_table);
    end
    
    % Channel Alignment
    % First make sure there's only 1 tile in each channel
    assert(max(path_table.x) == 1,...
        "More than 1 tile detected for marker %s. Images "+...
        "must be stitched first when multiple resolutions are present.",...
        config.markers(ncols == max(ncols)))
    assert(max(path_table.y) == 1,...
        "More than 1 tile detected for marker %s. Images "+...
        "must be stitched first when multiple resolutions are present.",...
        config.markers(nrows == max(nrows)))
    
    [config, path_table] = align_channels(config, path_table, equal_res);
end

%% Return outputs
if nargout == 1; path_table = []; end
if nargout == 0; config = []; end
    
end


function [config, path_table] = intensity_adjustment(config, path_table, nrows, ncols)
% Intensity adjustments

adj_file = fullfile(config.output_directory,"variables",'adj_params.mat');
loaded_params = false;

if ~isempty(config.update_intensity_channels)
    config.adjust_intensity(config.update_intensity_channels) = "update";
end

% Check for different selections between channels
if any(config.adjust_intensity == "true")
    % Check for previously calculated parameters
    if exist(adj_file,"file") == 2
        fprintf("%s\t Loading adjustment parameters \n",datetime('now'));
        load(fullfile(config.output_directory,'variables','adj_params.mat'),'adj_params')
        [adj_params, config] = check_adj_parameters(adj_params,config);
        
        % Additional checks on number of tiles
        for i = 1:length(adj_params.t_adj)
            if ~isempty(adj_params.t_adj{i})
            assert(size(adj_params.t_adj{i},1) == nrows(i), "Incorrect numner of rows "+...
                   "in the tile adjustment matrix for marker %s",config.markers(i))
            assert(size(adj_params.t_adj{i},2) == ncols(i), "Incorrect numner of columns "+...
                   "in the tile adjustment matrix for marker %s",config.markers(i))
            end
        end
        
        % Attach to config structure
        config.adj_params = adj_params;
    end
    
    % If all are loaded can return 
    if all(config.adjust_intensity == "true")
        return        
    else
        loaded_params = true;
    end
    
elseif all(config.adjust_intensity == "false")
    % No intensity adjustments
    fprintf("%s\t No intensity adjustments selected \n",datetime('now'));
    config.adj_params = [];

    if ~isempty(config.lowerThresh) && any(config.lowerThresh>1)
        idx = config.lowerThresh>1;
        config.lowerThresh(idx) = config.lowerThresh(idx)/65535;
    end

    if ~isempty(config.upperThresh) && any(config.upperThresh>1)
        idx = config.upperThresh>1;
        config.upperThresh(idx) = config.upperThresh(idx)/65535;
    end
    
    if ~isempty(config.signalThresh) && any(config.signalThresh>1)
        idx = config.signalThresh>1;
        config.signalThresh(idx) = config.signalThresh(idx)/65535;
    end
    return
    
elseif any(config.adjust_intensity ~= "false") && any(config.adjust_intensity ~= "update")
    error("Unrecognized selection for adjust_intensity. "+...
    "Valid options are ""true"", ""update"", or ""false"".")
end

if ~loaded_params || all(isequal(config.adjust_intensity,"update"))
    % Create new adj_params structure for all channels
    config.update_intensity_channels = 1:length(config.markers);

    % Define intensity adjustment parameters from scratch. Intensity
    % adjustment measurements should be made on raw images.
    fprintf("%s\t Defining new adjustment parameters \n",datetime('now'));        
    adj_fields = {'adjust_intensity','adjust_tile_shading',...
        'adjust_tile_position','lowerThresh','upperThresh','signalThresh'};

    % Update parameters from config structure
    for i = 1:length(adj_fields)
        adj_params.(adj_fields{i}) = config.(adj_fields{i});
    end
    
    % Save marker and img_location
    adj_params.markers = config.markers;
    if length(config.img_directory) == length(config.markers)
        adj_params.img_directory = config.img_directory;
    else
        adj_params.img_directory = repmat(config.img_directory,1,length(config.markers));
    end
end

% Calculate intensity adjustments
lowerThresh_measured = zeros(1,length(config.markers)); 
upperThresh_measured = lowerThresh_measured; signalThresh_measured = lowerThresh_measured; 
flatfield = cell(1,length(config.markers));  darkfield = flatfield;
for k = 1:length(config.markers)            
    if ~ismember(k,config.update_intensity_channels)
        continue
    end
    fprintf("%s\t Measuring Intensity for %s \n",datetime('now'),config.markers(k));
    stack = path_table(path_table.markers == config.markers(k),:);

    % Take measurements
    [lowerThresh_measured(k), upperThresh_measured(k), signalThresh_measured(k), t_adj, y_adj, flatfield{k}, darkfield{k}] = ...
        measure_images(config, stack, k);

    % Save new adjustments parameters structure
    if ~loaded_params
        adj_params.y_adj{k} = y_adj;
        adj_params.t_adj{k} = t_adj;
        adj_params.flatfield{k} = flatfield{k};
        adj_params.darkfield{k} = darkfield{k};
        adj_params.darkfield_intensity{k} = config.darkfield_intensity(k);
    end
end

% Check for user-defined intensity threshold values
[adj_params, config] = check_adj_parameters(adj_params, config,...
    lowerThresh_measured, upperThresh_measured, signalThresh_measured);
config.adj_params = adj_params;
config.adjust_intensity(config.adjust_intensity == "update") = "true";

% Save adjustment parameters to output directory
fprintf("%s\t Saving adjustment parameters \n", datetime('now'));
save(fullfile(config.output_directory,'variables','adj_params.mat'), 'adj_params')

% Save flatfield and darkfield as seperate variables
if isequal(config.adjust_tile_shading,"basic")
    fprintf("%s\t Saving flatfield and darkfield images \n",datetime('now'));
    save(fullfile(config.output_directory,'variables','flatfield.mat'), 'flatfield')
    save(fullfile(config.output_directory,'variables','darkfield.mat'), 'darkfield')
end

% Save samples
if isequal(config.save_samples,"true")
    fprintf('%s\t Saving samples \n',datetime('now'));    
    save_samples(config,'intensity',path_table)
end

end

function [config, path_table] = align_channels(config,path_table,equal_res)
% Channel Alignment

% Count number of x,y tiles
ncols = length(unique(path_table.x));
nrows = length(unique(path_table.y));
nb_tiles = ncols * nrows;
position_mat = reshape(1:nb_tiles,[ncols,nrows])';

% Check if z resolution is consistent for all channels
if ~equal_res
    error("NuMorph currently does not support multi-channel alignment at "+...
            "multiple z resolutions")
end

switch config.channel_alignment
    case 'translation'
        fprintf("%s\t Aligning channels by translation \n",datetime('now'));
        
        % Determine z displacement from reference channel
        % Check first if .mat file exists in output directory
        if exist(fullfile(config.output_directory,'variables','z_displacement_align.mat'),'file') == 2 && ...
                ~isequal(config.update_z_adjustment,"true")
            fprintf("%s\t Loading z displacement matrix \n",datetime('now'));
            load(fullfile(config.output_directory,'variables','z_displacement_align.mat'),'z_displacement_align')
        else
            for k = 2:length(config.markers)
                fprintf("%s\t Performing channel alignment z calculation for %s/%s \n",...
                    datetime('now'),config.markers{k},config.markers{1});
                z_tile = zeros(1,numel(position_mat));
                ave_signal = zeros(1,numel(position_mat));
                for idx = 1:numel(position_mat)
                    [y,x] = find(position_mat==idx);
                    path_ref = path_table(path_table.x==x & path_table.y==y & path_table.markers == config.markers(1),:);
                    path_mov = path_table(path_table.x==x & path_table.y==y & path_table.markers == config.markers(k),:);
                    % This measures the displacement in z for a given channel to the reference
                    [z_tile(idx),ave_signal(idx)] = z_align_channel(config,path_mov,path_ref,k);
                    fprintf("%s\t Predicted z displacement of %d for tile %d x %d \n",...
                        datetime('now'),z_tile(idx),y,x);
                end
                z_matrix = reshape(z_tile,[ncols, nrows])';
                z_displacement_align.(config.markers(k)) = z_matrix;
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
            alignment_table{y,x} = align_by_translation(config,path_align,z_displacement_align);
            save(save_path,'alignment_table')
            
        end
        fprintf("%s\t Alignment completed! \n",datetime('now'));
        
        % Change image directory to aligned directory so that subsequent
        % steps load these images
        if isequal(config.save_images,"true")
            config.img_directory = fullfile(config.output_directory,"aligned");
            path_table = path_to_table(config,"aligned");
        end
    case 'elastix'
        fprintf("%s\t Aligning channels using B-splines \n",datetime('now'));
        
        % Warn user if save_images flag is off. Can't apply params during
        % stitching
        if isequal(config.save_images,"false") && isequal(config.load_alignment_params,"false")
            warning("save_image flag is set to ""false"". Images need to be saved "+...
                "for downstream processing after elastix alignment.")
            pause(5)            
        end
        
        % Create variable
        alignment_params = cell(nrows,ncols);
        save_path = fullfile(config.output_directory,'variables','alignment_params.mat');
        if ~exist(save_path,'file') || isequal(config.load_alignment_params,"false")
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
        if isequal(config.load_alignment_params,"false") && isempty(config.signalThresh)
            % First check if adj_params exists in variables folder
            if exist(fullfile(config.var_directory,'align_thresholds.mat'),'file') == 2
                load(fullfile(config.var_directory,'align_thresholds.mat'),'thresholds')
            else
                fprintf("%s\t Lower and upper intensity thresholds are unspecified but are required "+...
                    "for accurate alignment. Measuring these now... \n",datetime('now'));
                for i = 1:length(config.markers)
                    [thresholds.lowerThresh(i), thresholds.upperThresh(i), thresholds.signalThresh(i)] = measure_images(config,path_table,i,true);
                end
                save(fullfile(config.var_directory,'align_thresholds.mat'),'thresholds')
            end
            % Save into config
            config.lowerThresh = thresholds.lowerThresh;
            config.singalThresh = thresholds.signalThresh;
            config.upperThresh = thresholds.upperThresh;
        end
        
        % If pre-aligning, check for and load calculated translation parameters
        if isequal(config.pre_align,"true")
            if exist(fullfile(config.var_directory,'alignment_table.mat'),'file') == 2
                load(fullfile(config.var_directory,'alignment_table.mat'),'alignment_table')
            else
                fprintf("%s\t Could not locate alignment table for pre-alignment \n",datetime('now'));
                alignment_table = [];
            end
        else
            alignment_table = [];
        end

        % Perform alignment
        for i = tiles_to_align
            [y,x] = find(position_mat==i);
            path_align = path_table(path_table.x==x & path_table.y==y,:);
            if ~isempty(alignment_table)
                config.alignment_table = alignment_table{y,x};
            end
            fprintf("%s\t Aligning channels to %s for tile 0%dx0%d \n",...
                datetime('now'),config.markers{1},y,x);
            
            alignment_params{y,x} = elastix_channel_alignment(config,path_align,true);
            save(save_path,'alignment_params','-v7.3')
            
            % Save samples
            if isequal(config.save_samples,"true")
                fprintf('%s\t Saving samples \n',datetime('now'));
                save_samples(config,'alignment',path_align)
            end
        end
                
        % Change image directory to aligned directory so that subsequent
        % steps load these images
        config.img_directory = fullfile(config.output_directory,'aligned');
        path_table = path_to_table(config,"aligned",false,true);
                
        % Update tile intensity adjustments using newly aligned images.
        % Also, set light sheet width adjustments + flatfield adjustments
        % to false as these were applied during the alignment step
        config.adjust_tile_shading = "false";
        config.adj_params.adjust_tile_shading = "false";
            
       % In case applying tile adjustments, re-calculate thresholds
        if isequal(config.adjust_intensity,"true") && isequal(config.adjust_tile_intensity,"true")
            for k = 1:length(config.markers)            
                fprintf("%s\t Updating intensity measurements for marker %s using "+...
                    "newly aligned images \n",datetime('now'),config.markers(k));
                stack = path_table(path_table.markers == config.markers(k),:);
                [~,~,t_adj] = measure_images(config, stack, k);
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

% Check number of loaded tiles
ncols = length(unique(path_table.x));
nrows = length(unique(path_table.y));
if nb_tiles ~= ncols*nrows
    warning("Number of aligned tiles and raw image tiles are not equal")
end

end


function [config, path_table] = perform_stitching(config,path_table)
% Stitching

% Count number of x,y tiles
ncols = length(unique(path_table.x));
nrows = length(unique(path_table.y));

if isequal(config.stitch_img,"true") || isequal(config.stitch_img, "load")
    fprintf("%s\t Perfoming iterative 2D image stitching \n",datetime('now'));

    % Check for aligned images
    if isequal(config.img_directory, fullfile(config.output_directory,'aligned'))
        fprintf("%s\t Using images from aligned directory \n",datetime('now'));
    elseif isequal(config.load_alignment_params,"true")
        fprintf("%s\t Loading channel alignment translations \n",datetime('now'));
        load(fullfile(config.output_directory,'variables','alignment_table.mat'),'alignment_table')
    end
    
    % Check for lowerThresh
    if isempty(config.lowerThresh)
        fprintf("%s\t Lower and upper intensity thresholds are unspecified but are required "+...
                "for stitching. Measuring these now... \n",datetime('now'));
        [config.lowerThresh, config.upperThresh] = measure_images(config,path_table, 1,true);
    end

    % Trim markers to only the ones specified
    if ~isempty(config.stitch_sub_channel)
       markers = config.markers(config.stitch_sub_channel);
       path_table = path_table(ismember(path_table.markers,markers),:);
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
        % Load z displacement info from a matrix
        fprintf("%s\t Loading z displacement matrix for stitching \n",datetime('now'));
        load(fullfile(config.var_directory, 'z_disp_matrix.mat'), 'z_disp_matrix')
        assert(nrows == size(z_disp_matrix,1) && ncols == size(z_disp_matrix,2), "Loaded z adjustment parameters do not match number of tiles "+...
            "detected in the input image directory. Check input image directory or update z adjustment")
        z_adj = apply_adjusted_z(path_table, z_disp_matrix);
        [~,z_idx] = setdiff(path_table.file,z_adj.file);
        path_table(z_idx,:) = []; 
        path_table.z_adj = z_adj.z_adj;
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
else
    fprintf("%s\t No image stitching selected \n",datetime('now'));
end

end
