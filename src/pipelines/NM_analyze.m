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
if isfield(config,'mri_directory') && ~isempty(config.mri_directory)
    [path_table, path_table_nii] = path_to_table(config);
else
    path_table = path_to_table(config);
    path_table_nii = [];
end

% If markers ignored, add these from raw
if ~isempty(path_table)
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
    if ntiles ~= 1
        warning("Check if use_processed_images field needs to be updated to stitched directory")
        assert(ntiles == 1, "To perform analysis, there should be only 1 tile for "+...
            "each channel in the image dataset.")
    end
end

%% Generate annotation .mat file if provided custom
%if isequal(config.use_annotation_mask,"true") && ~isempty(config.annotation_file)  
%    [~,b] = fileparts(config.annotation_file);
    %annot_file = fullfile(home_path,'data','masks',strcat(b,'.mat'));
%    annot_file = fullfile(config.output_directory,'variables',strcat(b,'.mat'));
    
%    if ~isfile(annot_file)
%        fprintf('%s\t Converting custom annotations \n',datetime('now'))
%        generate_annotations_from_file(config)
%    end
%end

%% Run single step and return if specified
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Note no specific checks on tiles/markers included
if nargin>1 && isequal(step,'resample') 
    config = perform_resampling(config, path_table);
    if nargout == 1; path_table = []; end
    if nargout < 1; clear config; end
    return
elseif nargin>1 && isequal(step,'register')
    config = perform_registration(config,path_table_nii);
    create_final_structure(config);
    if nargout == 1; path_table = []; end
    if nargout < 1; clear config; end
    return
elseif nargin>1 && isequal(step,'count')
    path_table = path_table(path_table.markers == config.markers(1),:);
    config = perform_counting(config,path_table);
    create_final_structure(config);
    if nargout == 1; path_table = []; end
    if nargout < 1; clear config; end
    return
elseif nargin>1 && isequal(step,'classify')
    config = perform_classification(config, path_table);
    create_final_structure(config);
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
    create_final_structure(config);
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
    if ~isempty(config.registration_channels)
        resample_channels = config.registration_channels;
    else
        resample_channels = unique(path_table.channel_num)';
    end
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
resample_path_table(path_table, config, resample_channels);

end


function config = perform_registration(config, path_table_mri)
% Image registration
home_path = fileparts(which('NM_config'));

if isequal(config.register_images,"false")
    % Check if I_mask exists
    mask_var = fullfile(config.output_directory,'variables',strcat(config.sample_id,'_mask.mat'));
    if ~isfile(mask_var) && ~isequal(config.use_annotation_mask,"false")
        error("Analysis is configured to use annotations however none exist. "+...
            "Run registration to generate an annotation mask or specify custom annotations")
    end
    return
end

I_mask = [];
if ~isempty(config.annotation_file)
    % Custom annotation
     if isequal(config.annotation_mapping,"image")
         % Matches images, save as .mat and return
         if isfile(config.annotation_file)
             I_mask = read_img(config.annotation_file);
             mask_var = fullfile(config.output_directory,...
                 'variables',strcat(config.sample_id,'_mask.mat'));
              fprintf('%s\t Saving custom annotations that are already mapped to images \n',datetime('now'))
             save(mask_var,'I_mask','-v7.3')
             return
         else
             error("Could not locate custom annotation file %s",config.custom_annotation_file)
         end
     else
         % Matches atlas file, still need to perform registration. Load
         % mask here
         I_mask = read_img(config.annotation_file);
     end
end
fprintf('%s\t Performing image registration \n',datetime('now'))

% Get which structures
[~,structures] = fileparts(config.structures);
r1_marker = config.markers(config.registration_channels);

% Get final direction with inverse
final_direction = config.registration_direction;
%if isequal(config.calculate_inverse,"true")
%    if isequal(final_direction,"atlas_to_image")
%        final_direction = "image_to_atlas";
%    elseif isequal(final_direction,"image_to_atlas")
%        final_direction = "atlas_to_image";
%    elseif isequal(final_direction,"mri_to_atlas")
%        final_direction = "altas_to_mri";
%    elseif isequal(final_direction,"atlas_to_mri")
%        final_direction = "mri_to_atlas";
%    elseif isequal(final_direction,"image_to_mri")
%        final_direction = "mri_to_image";
%    elseif isequal(final_direction,"mri_to_image")
%        final_direction = "image_to_mri";
%    end
%end

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

% Check atlas files
if ~loaded && contains(config.registration_direction,"atlas")
    atlas_path = arrayfun(@(s) {char(fullfile(config.home_path,'data','atlas',s))},config.atlas_file);
    assert(all(isfile(atlas_path)), "Could not locate Allen Reference Atlas .nii file specified")
    
    % Resize
    atlas_res = cellfun(@(s) repmat(str2double(regexp(s,'\d*','Match')),1,3),{config.atlas_file},'UniformOutput',false);

    % Rule: all atlas file must be at the same resolution
    assert(all(atlas_res{1} == atlas_res{end}),"All loaded atlas file must be at the same resolution")
end

% Perform pairwise registration
if ~loaded || isequal(config.register_images,"update")
    % Get moving image paths, subset channels, save into config
    if isequal(config.registration_direction,"image_to_atlas") ||...
            isequal(config.registration_direction,"image_to_mri")
        [~, resample_table] = path_to_table(config,'resampled');
        idx = ismember(resample_table.markers,[1,config.markers(config.registration_channels)]);
        mov_img_path = resample_table(idx,:).file;
        config.mov_res = resample_table.y_res;
        config.mov_orientation = config.orientation;
        config.mov_direction = "image";
        config.mov_channels = config.markers(config.registration_channels);
        if isequal(config.registration_prealignment,"image") ||...
                isequal(config.registration_prealignment,"both")
            config.mov_prealign = true;
        else
            config.mov_prealign = false;
        end
                
    elseif isequal(config.registration_direction,"mri_to_atlas") ||...
            isequal(config.registration_direction,"mri_to_image")
        mov_img_path = path_table_nii.file;
        config.mov_res = config.mri_resolution;
        config.mov_orientation = config.mri_orientation;
        config.mov_direction = "mri";
        config.mov_channels = config.mri_channels;
        if isequal(config.registration_prealignment,"mri") ||...
                isequal(config.registration_prealignment,"both")
            config.mov_prealign = true;
        else
            config.mov_prealign = false;
        end
        
    elseif isequal(config.registration_direction,"atlas_to_image") ||...
            isequal(config.registration_direction,"atlas_to_mri")
        % Type of atlas
        t = cell(1,length(config.atlas_file));
        for i = 1:length(config.atlas_file)
            t{i} = table(fullfile(home_path,'data','atlas',config.atlas_file),25,25,25,...
                'VariableNames',{'file','y_res','x_res','z_res'});
        end
        mov_img_path = {cat(1,t{:}).file};
        config.mov_res = 25;
        config.mov_orientation = "ail";
        config.mov_direction = "atlas";
        config.mov_channels = "atlas";
        config.mov_prealign = false;
        
    else
        error("Incorrect registration direction specified")
    end

    % Get reference image paths
    if isequal(config.registration_direction,"atlas_to_image") ||...
            isequal(config.registration_direction,"mri_to_image")
        [~, resample_table] = path_to_table(config,'resampled');
        % Subset channels to register
        idx = ismember(resample_table.markers,[1,config.markers(config.registration_channels)]);
        ref_img_path = resample_table(idx,:).file;
        config.ref_res = resample_table.y_res;
        config.ref_orientation = config.orientation;
        config.ref_direction = "image";
        config.ref_channels = config.markers(config.registration_channels);
        if isequal(config.registration_prealignment,"image") ||...
                isequal(config.registration_prealignment,"both")
            config.ref_prealign = true;
        else
            config.ref_prealign = false;
        end
        
    elseif isequal(config.registration_direction,"atlas_to_mri") ||...
            isequal(config.registration_direction,"image_to_mri")
        ref_img_path = path_table_nii.file;
        config.ref_res = config.mri_resolution;
        config.ref_orientation = config.mri_orientation;
        config.ref_direction = "mri";
        config.ref_channels = config.mri_channels;
        if isequal(config.registration_prealignment,"mri") ||...
                isequal(config.registration_prealignment,"both")
            config.ref_prealign = true;
        else
            config.ref_prealign = false;
        end
        
    elseif isequal(config.registration_direction,"image_to_atlas") ||...
            isequal(config.registration_direction,"mri_to_atlas")
        % Type of atlas
        t = cell(1,length(config.atlas_file));
        for i = 1:length(config.atlas_file)
            t{i} = table(fullfile(home_path,'data','atlas',config.atlas_file),25,25,25,...
                'VariableNames',{'file','y_res','x_res','z_res'});
        end
        ref_img_path = {cat(1,t{:}).file};
        config.ref_res = 25;
        config.ref_orientation = "ail";
        config.ref_direction = "atlas";
        config.ref_channels = "atlas";
        config.ref_prealign = false;
    end

    % Calculate registration parameters
    % Moving image to reference image
    reg_params = register_to_atlas(config, mov_img_path, ref_img_path);

    % Save registration parameters
    save(reg_file,'reg_params')
    fprintf('%s\t Registration completed! \n',datetime('now'))
end

% Save a copy of registered images
if isequal(config.save_registered_images,"true")
    fprintf('%s\t Transforming and saving registered images \n',datetime('now'))
    save_registered_images(config,reg_params)
end

% Now working on making an annotation mask
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
if ~isempty(I_mask)
    % Custom mask. Mask expected to match orientation of the registered
    % atlas
    
elseif isequal(config.use_annotation_mask,"true")
    % Generate mask for selected structures
    fprintf('%s\t Generating mask for selected structures \n',datetime('now'))
    I_mask = gen_mask(config.structures, config.hemisphere, config.orientation);
else
    fprintf('%s\t Not using annotations. Skipping mask generation. \n',datetime('now'))
    return
end

% Transform annotations
if ~isempty(reg_params) && ~isempty(reg_params.atlas_to_image)
    fprintf('%s\t Applying transformation to annotation mask \n',datetime('now'))

    % Adjust sizes and spacing
    I_mask = imresize3(I_mask,0.4,'Method','nearest');
    s = 1;    % Default: brain registration performed at 25um, mask at 10um
    size1 = reg_params.image_size;%reg_params.native_img_size;
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
% Detect cell/nuclei centroids

% Check for centroids structure
path_centroids = fullfile(config.var_directory,'centroids.mat');
if isequal(config.count_nuclei,"false")
    fprintf('%s\t No cell counting selected\n',datetime('now'))
    return
elseif isequal(config.count_nuclei,"true") && isfile(path_centroids)
    fprintf('%s\t Centroids already detected. Skipping cell counting\n',datetime('now'))
    return
end
   
% Append results to csv while progessing down stack
path_save = fullfile(config.output_directory, sprintf('%s_centroids.csv',config.sample_id));

if isequal(config.count_method,"hessian")
    % Detect nuclei centroids using hessian-based blob detector
    path_sub = path_table(path_table.markers == config.markers(1),:); 
    
    % Check for intensity thresholds
    if isempty(config.lowerThresh) || isempty(config.signalThresh) ||...
            isempty(config.upperThresh)
        markers = config.markers;
        config.markers = config.markers(1);
        config = check_for_thresholds(config,path_sub,true);
        config.markers = markers;
    end
    
    % Load annotation mask
    if isequal(config.use_annotation_mask,"true")
        config.mask_file = fullfile(config.output_directory,...
            'variables',strcat(config.sample_id,'_mask.mat'));
        load(config.mask_file,'I_mask')
    else
        I_mask = [];
    end
    
    % Run prediction
    predict_centroids_hessian(config,path_sub,path_save,I_mask);
    
elseif isequal(config.count_method,"3dunet")
    % Detect nuclei centroids using 3dunet (requires python + conda environment)
    config.path_save = path_save;
    config.img_directory = fileparts(path_table(path_table.channel_num ==1,:).file{1});
    predict_centroids_3dunet(config)
end

% Read centroids csv file and resave as MATLAB structure
centroids.coordinates = readmatrix(path_save);
save(path_centroids,'centroids')
delete(path_save)

end


function config = perform_classification(config, path_table)
% Classify Cell-Types  

% Check for centroids structure
path_classes = fullfile(config.var_directory,'classes.mat');
if isequal(config.classify_cells,"false")
    fprintf('%s\t No cell classification selected\n',datetime('now'))
    return
elseif isequal(config.classify_cells,"true") && isfile(path_classes)
    fprintf('%s\t Centroids already classified. Skipping cell classification\n',datetime('now'))
    return
end

% Subset markers to classify
if isempty(config.classify_channels)
    % Use only channels at the same resolution as nuclear channel
    c = 1:length(config.markers);
    config.classify_channels = c(cellfun(@(s) isequal(s,config.resolution{1}),config.resolution));
else
    
end
class_markers = config.markers(config.classify_channels);
path_table = path_table(ismember(path_table.markers,class_markers),:);

% Load centroids
path_centroids = fullfile(config.var_directory,'centroids.mat');
if ~isfile(path_centroids)
    error("Could not locate centroids file in the output directory");
else
    fprintf('%s\t Loading centroid list \n',datetime('now'))
    load(path_centroids,'centroids');
end

% Check if all channel intensity have been measured
if ~isfield(centroids,'intensities') || isequal(config.remeasure_centroids,'true') ||...
        size(centroids.intensities,2) ~= length(class_markers)
    % No intensities measured
    centroids.intensities = remeasure_centroids(centroids.coordinates,path_table,config,class_markers);
    centroids.intensities = round(centroids.intensities);
    save(path_centroids,'centroids')
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

% Name class file
path_classes = fullfile(config.var_directory,'classes.mat');
path_save = fullfile(config.output_directory, sprintf('%s_classes.csv',config.sample_id));

if isequal(config.classify_method,'gmm')
    % Not supported anymore. Might add later

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
    
    % Create new patches for manual labeling or load a previous set
    if ~isfile(itable) || ~isfile(ftable)
        fprintf('%s\t Generating new centroid patches for training SVM model \n',...
            datetime('now'))
        get_centroid_patches(centroids,path_table,config);
    elseif isequal(config.load_patches,"true")
        get_centroid_patches('load',path_table,config);
    end
    
    % Get full patch features
    if ~isfile(ftable_full)
        fprintf("%s\t Generating complete table of centroid patch features. "+....
            "This may take some time...\n",datetime('now'));
        stable = get_patch_features(centroids,path_table,config);
        writematrix(stable,ftable_full)
    else
        fprintf('%s\t Loading complete table of patch features \n',datetime('now'))
        stable = readmatrix(ftable_full);
        if size(stable) ~= size(centroids,1)
           error("Number of patch feature table does not match the number of centroids loaded")
        end
    end
    
    % Check for all tables prior to running classification
    if isfile(itable) && isfile(ftable) && ~isfile(ctable)
        error("Could not locate annotated classifications in classifier directory. "+...
            "Generate annotations using generate_classifications(config)")
    else
        ftable = readtable(ftable);
        ctable = readtable(ctable);
        
        % Check size
        if size(ftable,1) ~= height(ctable)
            error("Size of annotated classifications do not match size of feature table. "+...
                "Re-annotate patches using current settings.")
        end
        
        % Check if any empty
        if any(ctable{:,1} == 0)
            error("Loaded annotations contain 0 values. All patches must be annotated with a value 1-9.")
        end
    end
    
    % Classify cells
    ct = classify_cells_svm(centroids, stable, config);
    
elseif isequal(config.classify_method,'threshold')
    % Use simple thresholding for predictions
    ct = classify_cells_threshold(centroids, config);

else
    error('%s\t Invalid classification method specified \n',datetime('now'))
end

% Construct new csv file for cell classes
cen_classes = cat(2,centroids.coordinates,centroids.annotations,ct);

% Remove outliers
if isempty(config.keep_classes)
    config.keep_classes = unique(ct(ct>0)');
end
keep_idx = ismember(ct,config.keep_classes);
cen_classes = cen_classes(keep_idx,:);
centroids.intensities = centroids.intensities(keep_idx,:);

% Print results
fprintf('\t Total centroids retained: %d\n', size(cen_classes,1))
u_classes = unique(cen_classes(:,5));
for i = 1:length(u_classes)
    idx = cen_classes(:,5) == u_classes(i);
    count = sum(idx);
    intensities = round(mean(centroids.intensities(idx,:),1));
    fprintf('\t Classified %d cells or %.3f of centroids as type %d.\n',...
        count,count/size(cen_classes,1),u_classes(i))
    fprintf('\t Average intensities: %s.\n',num2str(intensities))
end

% Save to structure
classes.centroids = cen_classes(:,1:3);
classes.annotations = cen_classes(:,4);
classes.classes = cen_classes(:,5);
save(path_classes,'classes')

fprintf('%s\t Classification of sample %s completed! \n',...
        datetime('now'),config.sample_id)
end

