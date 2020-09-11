function NM_analyze(config)
%--------------------------------------------------------------------------
% NuMorph analysis pipeline designed to perform image registration, nuclei
% counting, and cell-type classification based on nuclear protein markers
% in whole brain images.
%--------------------------------------------------------------------------

if nargin<1 
    load 'NM_variables.mat'
    config = load(fullfile('templates','NM_variables.mat'));
end
fprintf('%s\t Working on sample %s \n',datetime('now'),config.sample_name)

%% Read Filename Information
% Load image paths for registration
if any(strcmp([resample_image,register_image],"true"))
    if isequal(resample_image,"true")
        % Resample images to desired atlas resolution
        fprintf('%s\t Reading image filename information from stitched directory \n',datetime('now'))
        path_cell{1} = dir(img_directory);
        if isequal(fullfile(output_directory,'stitched'),img_directory)
            location = 'stitched';
            path_table_stitched = path_to_table(path_cell,location,markers,channel_num);
        end

        % Unless specified otherwise, resample all markers
        if isempty(config.resample_markers)
            resample_markers = 1:length(markers);
        end
        for k = resample_markers
            fprintf('%s\t Resampling channel %s\n',datetime('now'),markers(k))            
            path_sub = path_table_stitched(path_table_stitched.markers == markers(k),:);
            resample_img_to_atlas(path_sub, output_directory, resolution, resample_res);
        end
        resample_image = "load";
    end

    if exist(fullfile(output_directory,'resampled'),'dir') == 7 && isequal(resample_image,"load")
        % Load images from resampled directory
        fprintf('%s\t Reading image filename information from resampled directory \n',datetime('now'))
        path_cell{1} = dir(fullfile(output_directory,'resampled'));
        location = 'resampled';
        path_table_resampled = path_to_table(path_cell,location,markers,channel_num);        
    else
        error('%s\t Could not locate resampled directory in the specified image directory. ' +...
            "Set resample_image to ""true"" to generate images for registration",string(datetime('now')));
    end
end

%% Function for Registration
switch register_image
    case 'true'
        mov_img_path = path_table_resampled.file{1};
        
        % Calculate registration parameters
        reg_params = register_to_atlas(config, mov_img_path);

        % Apply transformation to other channels in the sample
        if isequal(config.save_registered_image,'true')
            for i = 2:height(path_table_resampled)
                fprintf('%s\t Applying transform to %s images\n',datetime('now'),markers(i)); 
                mov_img_path = path_table_resampled.file{i};
                apply_transform_to_resampled(mov_img_path,reg_params)
            end
        end
    
    case 'load'
        % Attempt to load registration parameters
        fprintf('%s\t Loading registration parameters \n',datetime('now'))
        try
            load(fullfile(output_directory,'variables','reg_params.mat'),'reg_params')
        catch ME
            error('%s\t Could not locate registration parameters \n',string(datetime('now')))
        end
    case 'false'
        fprintf('%s\t No image registration selected \n',datetime('now'))
        reg_params = [];
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
        
        measure_dice = true;
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
        path_table_stitched = path_to_table(path_cell,location,[],[]);
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

%% Generate Background Mask 
if isequal(measure_local_background,"true")
    fprintf('%s\t Generating image of local background intensities \n',datetime('now'))
    [B, S] = create_background_image(path_table_stitched, config);
elseif isequal(measure_local_background,"load")
    fprintf('%s\t Loading previous background images \n',datetime('now'))
    path_background = fullfile(output_directory, 'background');
    files = dir(path_background);
    B = cell(1,length(markers));
    S = cell(1,length(markers));
    for i = 1:length(markers)
        img_idx = contains({files.name}, "background.nii") &...
            contains({files.name}, markers(i));
        if any(img_idx)
            B{i} = niftiread(fullfile(path_background,files(img_idx).name));
        end
        img_idx = contains({files.name}, "back_std.nii") &...
                       contains({files.name}, markers(i)); 
        if any(img_idx)
           S{i} = niftiread(fullfile(path_background,files(img_idx).name));
        end
    end
else
    B = cell(1,length(markers));
    S = cell(1,length(markers));
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
    centroids = remeasure_centroids(centroids,path_table_stitched,config);
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
        ct = classify_cells_svm(centroids, path_table_stitched, config);
        
    case 'threshold'
        ct = classify_cells_threshold(centroids, config, B, S);
        centroids(:,end+1) = ct;
    case 'none'
        fprintf('%s\t No centroid list so skipping cell classification \n',datetime('now'))
end

fprintf('%s\t Analysis steps completed! \n',datetime('now'))

end