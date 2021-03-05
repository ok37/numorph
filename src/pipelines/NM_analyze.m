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
home_path = fileparts(which('NM_config'));

% Default to run full pipeline
if nargin<2
    step = 'analyze';
end

% Make an output directory
if ~isfolder(config.output_directory)
    mkdir(config.output_directory);
end

% Make a variables directory
if ~isfolder(config.var_directory)
    mkdir(config.var_directory);
end

fprintf("%s\t Working on sample %s \n",datetime('now'),config.sample_id)

%% Read image filename information
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate table containing image information
path_table = path_to_table(config);

% If markers ignored, add these from raw
if ~all(ismember(config.markers,unique(path_table.markers)))
    idx = find(~ismember(config.markers,unique(path_table.markers)));
    
    for i = 1:length(idx)
        config2 = config;
        config2.markers = config.markers(idx(i));
        config2.channel_num = config.channel_num(idx(i));
        try 
            path_table = vertcat(path_table,path_to_table(config2,'raw',false));
            path_table.channel_num = arrayfun(@(s) find(s == config.markers),path_table.markers);
        catch
            warning("Could not load marker %s ignored from processing.",config.markers(idx(i)));
        end
    end
    clear config2
end

% Count number of x,y tiles for each channel
ntiles = length(unique([path_table.x,path_table.y]));
assert(ntiles == 1, "To perform analysis, there should be only 1 tile for "+...
    "each channel in the image dataset.")

%% Generate annotation .mat file if provided custom
if ~ismember(config.use_annotation_mask,["true","false"])    
    [~,b] = fileparts(config.use_annotation_mask);
    annot_file = fullfile(home_path,'data','masks',strcat(b,'.mat'));

    fprintf('%s\t Using custom annotations \n',datetime('now'))
    if ~isfile(annot_file)
        generate_annotations_from_file(config)
    end
end

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
    path_table = path_table(path_table.markers == config.markers(1),:);
    config = perform_counting(config,path_table);
    if nargout == 1; path_table = []; end
    if nargout < 1; clear config; end
    return
elseif nargin>1 && isequal(step,'classify')
    config = perform_classification(config, path_table);
    if nargout == 1; path_table = []; end
    if nargout < 1; clear config; end
    return
end

%% Run full analysis pipeline
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isequal(step,'analyze') 
    config = perform_resampling(config, path_table);
    config = perform_registration(config);
    config = perform_counting(config,path_table);
    config = perform_classification(config, path_table);
    fprintf('%s\t Analysis steps completed! \n',datetime('now'))
end

end


function config = perform_resampling(config, path_table)
% Image resampling

res = config.resample_resolution;
res_path = fullfile(config.output_directory,'resampled');
config.resampled_paths = cell(1,length(config.markers));
flags = true(1,length(config.markers));

% Get channels to resample
if isempty(config.resample_channels)
    resample_channels = unique(path_table.channel_num)';
    flags(~resample_channels) = false;
else
    resample_channels = config.resample_channels;
end

% Create resampled directory
if ~isfolder(res_path)
    mkdir(res_path)
else
    if ~isequal(config.resample_images,"update")
        % Check for images alread resampled
        for i = resample_channels
            filename = fullfile(res_path,sprintf('%s_C%d_%s_%d.nii',...
                config.sample_id,i,config.markers(i),res));
            if ~isfile(filename)
                flags(i) = false;
            end
        end
        
        % If all channels resampled at correct resolution, return
        if ~isequal(config.resample_images,"update")
            if all(flags)
                fprintf('%s\t All selected channels have been resampled\n',datetime('now'))            
                return
            end
        end
    end
end

% Perform resampling
resample_path_table(path_table, config,resample_channels);

end


function config = perform_registration(config)
% Image registration

if isequal(config.register_images,"false")
    % Check if I_mask exists
    mask_var = fullfile(config.output_directory,'variables',strcat(config.sample_id,'_mask.mat'));
    if ~isfile(mask_var) && ~isequal(config.use_annotation_mask,"false")
        error("Analysis is configured to use annotations however none exist. "+...
            "Run registration to generate an annotation mask or specify custom annotations")
    end
    return
end
fprintf('%s\t Performing image registration \n',datetime('now'))

% Get which structures
[~,structures] = fileparts(config.structures);
r1_marker = config.markers(config.register_channels);
home_path = fileparts(which('NM_config'));

% Get final direction with inverse
final_direction = config.direction;
if isequal(config.calculate_inverse,"true")
    if isequal(final_direction,"atlas_to_image")
        final_direction = "image_to_atlas";
    elseif isequal(final_direction,"image_to_atlas")
        final_direction = "atlas_to_image";
    elseif isequal(final_direction,"mri_to_atlas")
        final_direction = "altas_to_mri";
    elseif isequal(final_direction,"atlas_to_mri")
        final_direction = "mri_to_atlas";
    elseif isequal(final_direction,"image_to_mri")
        final_direction = "mri_to_image";
    elseif isequal(final_direction,"mri_to_image")
        final_direction = "image_to_mri";
    end
end

% Attempt to load registration parameters
reg_params = [];
loaded = false;
reg_file = fullfile(config.output_directory,'variables','reg_params.mat');
if isfile(reg_file) 
    % Load previous parameters if not updating
    if ~isequal(config.register_images,"update")
        fprintf('%s\t Loading previosuly calculated registration parameters \n',datetime('now'))
        load(reg_file,'reg_params')
    end
    
    % Check if loaded registration parameters contain the direction
    % specified in config
    if isfield(reg_params,final_direction)
        loaded = true;
    end
    if ~isequal(config.register_images,"update")
        fprintf('%s\t Parameters already exist for specified direction. Skipping registration\n',...
            datetime('now'))
    end
end

% Read file and perform registration
if ~loaded || isequal(config.register_images,"update")
    % Get resampled paths
    [~,resample_table] = path_to_table(config,'resampled');

    % Subset channels to register
    idx = ismember(resample_table.markers,[1,config.markers(config.register_channels)]);
    resample_table = resample_table(idx,:);

    % Calculate registration parameters
    reg_params = register_to_atlas(config, resample_table.file);

    % Save registration parameters
    save(reg_file,'reg_params')
    
    fprintf('%s\t Registration completed! \n',datetime('now'))
end

% Save a copy of registered images
reg_dir = fullfile(config.output_directory,'registered');
if isequal(config.save_registered_images,"true")
    fprintf('%s\t Transforming and saving registered images \n',datetime('now'))
    if ~isfolder(reg_dir)
        mkdir(reg_dir)
    end
    
    % Get resampled paths
    [~,resample_table] = path_to_table(config,'resampled');

    % Subset channels to register
    idx = ismember(resample_table.markers,config.markers(config.register_channels));
    resample_table = resample_table(idx,:);
        
    atlas_path = fullfile(home_path,'data','atlas',config.atlas_file(1));
    final_dest = strsplit(final_direction,"_");
    
    % Move final target image to registered dir
    if isequal(final_dest(3),"atlas")
        [~,fname] = fileparts(atlas_path);
        copyfile(atlas_path,strcat(reg_dir,'_target.nii'))
    elseif isequal(final_dest(3),"image")
        [~,fname] = fileparts(resample_table.file{1});
        copyfile(resample_table.file{1},fullfile(reg_dir,strcat(fname,'_target.nii')))
    end
    
    % Transform image to be registered
    if isequal(final_dest(1),"atlas")
        I_reg = niftiread(atlas_path);
        I_reg = transformix(I_reg,reg_params.atlas_to_image,[1,1,1], []);
        [~,fname] = fileparts(atlas_path);
    elseif isequal(final_dest(1),"image")
        I_reg = niftiread(resample_table.file{1});
        I_reg = transformix(I_reg,reg_params.image_to_atlas,[1,1,1], []);
        [~,fname] = fileparts(resample_table.file{1});
    end
    fname = fullfile(reg_dir,strcat(fname,'_registered.nii'));
    niftiwrite(uint16(I_reg),fname)
end

% Name mask image based on direction
if isequal(final_direction,"atlas_to_image") || isequal(final_direction,"mri_to_image")
    annot_marker = r1_marker;
elseif isequal(final_direction,"image_to_atlas") || isequal(final_direction,"image_to_atlas")
    [~,annot_marker] = fileparts(config.atlas_file);
elseif isequal(final_direction,"image_to_mri") || isequal(final_direction,"atlas_to_mri")
    annot_marker = config.mri_markers(1);
end

annot_file = fullfile(reg_dir,sprintf('%s_%s_%d_%s.nii',...
            config.sample_id,strjoin(annot_marker,"_"),....
            config.resample_resolution,structures));

% Permute mask to match atlas file
if isequal(config.use_annotation_mask,"true")
    % Generate mask for selected structures
    fprintf('%s\t Generating mask for selected structures \n',datetime('now'))
    I_mask = gen_mask(config.hemisphere, config.structures);
    I_mask = permute_orientation(I_mask,'sla','ail');
elseif ~isequal(config.use_annotation_mask,"false")
   %I_mask = read_img(); 
    
end

% Generate annotations
if ~isempty(reg_params) && ~isempty(reg_params.atlas_to_image)
    fprintf('%s\t Applying transformation to annotation mask \n',datetime('now'))

    % Adjust sizes and spacing
    I_mask = imresize3(I_mask,0.4,'Method','nearest');
    s = 1;    % Default: brain registration performed at 25um, mask at 10um
    size1 = reg_params.native_img_size;%reg_params.native_img_size;
    for j = 1:length(reg_params.atlas_to_image.TransformParameters)
        reg_params.atlas_to_image.TransformParameters{j}.FinalBSplineInterpolationOrder = 0;
        reg_params.atlas_to_image.TransformParameters{j}.Size = size1;
        reg_params.atlas_to_image.TransformParameters{j}.Spacing = [s,s,s];
    end

    % Transform atlas
    I_mask = transformix(I_mask,reg_params.atlas_to_image,[s,s,s], []);
    
    % Save mask as variable
    save_name = fullfile(config.output_directory,'variables',strcat(config.sample_id,'_mask.mat'));
    save(save_name,'I_mask','-v7.3')

    % Save a copy in registered directory
    if isequal(config.save_registered_images,"true")
        fprintf('%s\t Saving annotation mask to registered directory \n',datetime('now'))
        if ~isfolder(reg_dir)
            mkdir(reg_dir)
        end
        I_mask_small = imresize3(I_mask,size(I_reg),...
            'Method','nearest');
        niftiwrite(uint16(I_mask_small),annot_file)
    end
    
    fprintf('%s\t Annotations generated! \n',datetime('now'))
else
    error('%s\t No registration parameters to transform annotation mask',string(datetime('now')))
end

% Measure structure volumes
fprintf('%s\t Measuring structure volumes \n',datetime('now'))
measure_structure_volumes(config);
        
end


function config = perform_counting(config,path_table)
% Segment and/or Count Cells

path_centroids = fullfile(config.output_directory, sprintf('%s_centroids1.csv',...
    config.sample_id));

% Update centroid index if file already exists
if isfile(path_centroids)
    [path,fname,ext] = fileparts(path_centroids);
    fname = char(fname);
    fname(end) = num2str(str2double(fname(end))+1);
    path_save = fullfile(path,strcat(fname,ext));
else
    path_save = path_centroids;
end

if isequal(config.count_method,"hessian")
    % Detect nuclei centroids using hessian-based blob detector
    path_sub = path_table(path_table.markers == config.markers(1),:); 
    
    % Load annotation mask
    if isequal(config.use_annotation_mask,"true")
        config.mask_file = fullfile(config.output_directory,...
            'variables',strcat(config.sample_id,'_mask.mat'));
        load(config.mask_file,'I_mask')
    else
        I_mask = [];
    end
    
    % Run prediction
    predict_centroids_hessian(config,path_sub,I_mask);
    
elseif isequal(config.count_method,"3dunet")
    % Detect nuclei centroids using 3dunet (requires python + conda environment)
    config.path_save = path_save;
    config.img_directory = fileparts(path_table(path_table.channel_num ==1,:).file{1});
    predict_centroids_3dunet(config)
    
end


end


function config = perform_classification(config, path_table)
% Classify Cell-Types  

markers = config.markers(config.classify_channels);
path_table = path_table(ismember(path_table.markers,markers),:);

% Load previous cell counting results. Should be located in the output directory.
path_centroids = dir(config.output_directory);
path_centroids = path_centroids(arrayfun(@(s) startsWith(s.name,sprintf("%s_centroids",config.sample_id)) ...
    && endsWith(s.name,".csv"),path_centroids));

% Take most recent one if multiple
if length(path_centroids)>1
    s1 = arrayfun(@(s) strsplit(s.name,{'centroids','.csv'}),path_centroids,...
        'UniformOutput',false);
    s1 = cellfun(@(s) str2double(s(2)),s1);
    path_centroids = path_centroids(s1 == max(s1));
end

% Load centroids
if isempty(path_centroids)
    error("Could not locate centroids file in the output directory");
else
    fprintf('%s\t Loading centroid list from %s \n',datetime('now'),path_centroids.name)
    path_centroids = fullfile(path_centroids.folder,path_centroids.name);
    centroids = readmatrix(path_centroids);
end

% Check if all channel intensity have been measured
ncols = 4 + length(markers);

if size(centroids,2) == 4 || isequal(config.remeasure_centroids,'true')
    % No intensities measured
    m_int = remeasure_centroids(centroids,path_table,config,markers);    
    centroids(:,5:ncols) = m_int;
    writematrix(centroids,path_centroids);
    
elseif size(centroids,2) == 5 && length(markers)>1
    % Reference channel already measured
    m_int = remeasure_centroids(centroids,path_table,config,markers(2:end));
    centroids(:,6:ncols) = m_int;
    writematrix(centroids,path_centroids);

elseif size(centroids,2) ~= ncols
    % Incorrect number of columns in centroid file
    error("Error in centroids .csv file. Number of data columns does not match"+...
        " expected number of columns based on number of markers specified")
end

% Calculate automatic minimum threshold if unspecified
if isempty(config.min_intensity)
    if ~isfile(fullfile(config.output_directory,'variables','thresholds.mat'))
        [~,~,thresh] = measure_images(config, path_table, 1, true);
    else
        thresh = load_var(config,'thresholds');
        thresh = thresh.signalThresh(1);
    end
    config.min_intensity = thresh;
    fprintf('%s\t Using minimum centroid intensity: %.0f \n',datetime('now'),thresh*65535)
end

% Subset channels from path_table to classify
if isempty(config.classify_channels)
    c = 1:length(markers);
    c = c(cellfun(@(s) isequal(s,config.resolution{1}),config.resolution));
else
    c = config.classify_channels;
end
path_table = path_table(ismember(path_table.channel_num,c),:);

% Name class file
path_classes = fullfile(config.output_directory, sprintf('%s_classes1.csv',...
    config.sample_id));
if isfile(path_classes)
    [path,fname,ext] = fileparts(path_classes);
    fname = char(fname);
    fname(end) = num2str(str2double(fname(end))+1);
    path_save = fullfile(path,strcat(fname,ext));
else
    path_save = path_classes;
end

if isequal(config.classify_method,'gmm')
    [ct, p, gm] = classify_cells_gmm(centroids, config);
    centroids(:,8) = ct;

    %path_centroids = fullfile(output_directory, sprintf('%s_centroids.csv',sample_id));
    %writematrix(centroids,path_centroids);
        
    %I_final = create_image_slices(centroids, path_table_stitched, config);
    %I_final = create_image_slice(centroids, path_table_stitched, config);
        
    %if isequal(config.save_counts,'true') || isequal(config.save_counts,'overwrite')
    %    fprintf('%s\t Saving centroid list \n',datetime('now'))
    %    path_centroids = fullfile(output_directory, sprintf('%s_centroids2_new.csv',sample_id));
    %    writematrix(centroids_new,path_centroids);
    %end
elseif isequal(config.classify_method,'svm')
    % Classify cell-types by trained SVM
    
    % Check for centroid patches    
    save_directory = fullfile(config.output_directory,'classifier');
    if ~isfolder(save_directory)
        mkdir(save_directory)
    end
    
    % Read classifications + patch feautre info
    itable = fullfile(config.output_directory,'classifier',...
        sprintf('%s_patch_info.csv',config.sample_id));
    ftable = fullfile(config.output_directory,'classifier',...
        sprintf('%s_patch_features.csv',config.sample_id));
    ctable = fullfile(config.output_directory,'classifier',...
        sprintf('%s_classifications.csv',config.sample_id));
    ftable_full = fullfile(config.output_directory,'classifier',...
        sprintf('%s_full_patch_features.csv',config.sample_id));
    
    if ~isfile(itable) || ~isfile(ftable)
        fprintf('%s\t Generating new centroid patches for training SVM model \n',...
            datetime('now'))
        get_centroid_patches(centroids,path_table,config);
    elseif isequal(config.load_patches,"true")
        %get_centroid_patches('load',path_table,config);
    end
    if isfile(itable) && isfile(ftable) && ~isfile(ctable)
        %error("Could not locate annotated classifications in classifier directory."+...
        %    "Generate annotations using generate_annotation(config)")
    end
    
    % Get full patch features
    if ~isfile(ftable_full)
        fprintf('%s\t Generating complete table of centroid patch features \n',...
            datetime('now'))
        stable = get_patch_features(centroids,path_table,config);
        writematrix(stable,ftable_full)
        return
    else
        fprintf('%s\t Loading complete table of patch features \n',datetime('now'))
        stable = readmatrix(ftable_full);
        if size(stable) ~= size(centroids,1)
           error("Number of patch feature table does not match the number of centroids loaded")
        end
    end
    
    % Classify cells
    ct = classify_cells_svm(centroids, stable, config);
    
    % Construct new csv file for cell classes
    cen_classes = centroids(:,1:4);
    cen_classes = horzcat(cen_classes,ct);
    
    % Remove outlers
    r_idx = ismember(ct,config.keep_classes);
    cen_classes = cen_classes(r_idx,:);
    
    % Write matrix
    writematrix(cen_classes,path_save)
    
elseif isequal(config.classify_method,'threshold')
    ct = classify_cells_threshold(centroids, config, B, S);
    centroids(:,end+1) = ct;
else
    error('%s\t Invalid classification method specified \n',datetime('now'))
end

% Print results
fprintf('\t Total centroids retained: %d\n', size(cen_classes,1))
classes = unique(cen_classes(:,5));
for i = 1:length(classes)
    count = sum(cen_classes(:,5) == classes(i));
    fprintf('\t Classified %d cells or %.3f of centroids as type %d.\n',...
        count,count/size(cen_classes,1),classes(i))
end

fprintf('%s\t Classification of sample %s completed! \n',...
        datetime('now'),config.sample_id)
end