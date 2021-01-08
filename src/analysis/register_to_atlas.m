function reg_params = register_to_atlas(config, mov_img_path, num_points, mov_mask, atlas_mask)
% Register images to the reference atlas using elastix via melastix
% wrapper. Note: the final registration parameters are stored in the
% reg_params variable. While it is simpler to just register the atlas to
% the image, you may, in some cases, get better accuracy by registering the
% image to the atlas and calculating the inverse parameters. 

% Defaults
high_prct = 99;
Gamma = 0.9;

% Unpack config variables
resample_res = config.resample_resolution;
registration_method = config.registration_method;
hemisphere = config.hemisphere;
orientation = config.orientation;
output_directory = config.output_directory;
mask_coordinates = config.mask_coordinates;
direction = config.direction;
calc_inverse = config.calculate_inverse;
points_file = config.points_file;
save_registered_image = config.save_registered_image;
atlas_file = config.atlas_file;
perm = [0,0];

% Unless testing, use all points
if nargin == 2
    num_points = 'all';
end

% Convert to cell
if ~iscell(mov_img_path)
    mov_img_path = {mov_img_path};
end

% Mung atlas file
if isfile(atlas_file)
    atlas_path = {atlas_file};
else
    atlas_path = fullfile(config.home_path,'data','atlas',atlas_file); 
    if ~isfile(atlas_path)
        error("Could not locate Allen Reference Atlas .nii file specified")
    end
end

% Load atlas + moving image 
fprintf('%s\t Reading resampled and atlas images\n',datetime('now')); 

% MATLAB will usually rotate and flip axis dimensions based on the data
% type in the niftifile
mov_img_final = cell(1,length(mov_img_path));
for i = 1:length(mov_img_path)
    mov_img = niftiread(mov_img_path{i});
    if isequal(class(mov_img),'double')
        mov_img = imrotate(mov_img,90);
        mov_img = flip(mov_img,1);
    else
       mov_img = double(mov_img); 
    end
    mov_img_final{i} = mov_img;
end

% Load ARA
atlas_img_final = cell(1,length(atlas_path));
for i = 1:length(atlas_path)
    atlas_img = niftiread(atlas_path{i});
    % Transform atlas_img based on which sample is being imaged
    if isequal(hemisphere, "right")
        perm(1) = 1;
        atlas_img = flip(atlas_img,1);
    elseif isequal(hemisphere,"whole")
        perm(1) = 2;
        atlas_img2 = flip(atlas_img,3);
        atlas_img = cat(3,atlas_img,atlas_img2);
        if isequal(orientation,"lateral")
            warning('%s\t Using lateral orientation on whole brain image\n',datetime('now'));
        end
    end
    if isequal(orientation,"dorsal")
        perm(2) = 1;
        atlas_img = permute(atlas_img,[1,3,2]);
        atlas_img = flip(atlas_img,3);
    elseif isequal(orientation,"ventral")
        perm(2) = 2;
        atlas_img = permute(atlas_img,[1,3,2]);
    end
    atlas_img_final{i} = atlas_img;
end

% Resize
native_img_size = size(mov_img);
atlas_res = repmat(str2double(regexp(atlas_file,'\d*','Match')),1,3);
res_adj = atlas_res./resample_res;
if res_adj(1) > 1
    % Atlas is at a lower resolution. Downsample image
    mov_img_final = cellfun(@(s) imresize3(s,res_adj.*size(s)),mov_img_final);
elseif res_adj(1) < 1
    % Image is at a lower resolution. Downsample atlas
    atlas_img_final = cellfun(@(s) imresize3(s,res_adj.*size(s)),atlas_img_final);
end

% Rescale intensities
mov_img_final = cellfun(@(s) double(imadjustn(uint16(s),...
    [0,prctile(s(:),high_prct)/65535],[],Gamma)), mov_img_final,...
    'UniformOutput', false);

% Add mask if coordinates are specified
mask = [];
if ~isempty(mask_coordinates)
    mask = zeros(size(mov_img));
    mask(mask_coordinates(1):mask_coordinates(2),:,:) = 1;
end
    
% Chain registration parameters. Parameter files are found in
% ./elastix_parameter_files/atlas_registration. Update specific parameters
% by editting these files.
parameter_location = fullfile(config.home_path,'elastix_parameter_files','atlas_registration');
parameter_paths = cell(1,length(registration_method));
use_points = false;
for i = 1:length(registration_method)
    reg_type = registration_method(i);
    switch reg_type
        case 'a'
            fprintf('%s\t Performing affine registration\n',datetime('now'));
            parameter_paths{i} = fullfile(parameter_location,'ElastixParameterAffine.txt');
        case 'b'
            fprintf('%s\t Performing b-spline registration\n',datetime('now')); 
            parameter_paths{i} = fullfile(parameter_location,'ElastixParameterBSpline.txt');

        case 'p'
            fprintf('%s\t Performing points registration\n',datetime('now'));
            
            % Load points from BigWarpJ .csv file
            [mov_points,atlas_points] = load_points_from_bdv(output_directory, points_file);
            
            % Trim points if not using all
            if isequal(num_points,'all')
                num_points = size(mov_points,1);
            end
            mov_points = mov_points(1:num_points,:);
            atlas_points = atlas_points(1:num_points,:);
            use_points = true;                        
            
            parameter_paths{1} = fullfile(parameter_location,'ElastixParameterAffinePoints.txt');
            parameter_paths{2} = fullfile(parameter_location,'ElastixParameterPoints.txt');
            break
    end
end

% Create temporary directory for saving images
outputDir = fullfile(config.output_directory,sprintf('tmp_reg_%d',randi(1E4)));
if ~exist(outputDir,'dir')
    mkdir(outputDir)
end

reg_params = struct('img_to_atlas',[],'atlas_to_img',[]);

% Perform registration
if isequal(direction,'atlas_to_image')
    if isempty(mask)
        % Register atlas to image with or without corresponding points
        if ~use_points
            [reg_params.atlas_to_img,reg_img]=elastix(atlas_img,mov_img,...
                outputDir,parameter_paths,'threads',[]);
        else
            [reg_params.atlas_to_img,reg_img]=elastix(atlas_img,mov_img,...
                outputDir,parameter_paths,'fp',mov_points,'mp',atlas_points,'threads',[]);
        end
    else
        % Register atlas to image with or without corresponding points +
        % with a mask
        if ~use_points
            [reg_params.atlas_to_img,reg_img]=elastix(atlas_img,mov_img,...
                outputDir,parameter_paths,'mMask',mov_mask,'fMask',atlas_mask,'threads',[]);
        else
            [reg_params.atlas_to_img,reg_img]=elastix(atlas_img,mov_img,...
                outputDir,parameter_paths,'fp',mov_points,'mp',atlas_points,...
                'mMask',mov_mask,'fMask',atlas_mask,'threads',[]);
        end
        % Remove transformed image from structure
        reg_params.atlas_to_img.transformedImages = [];
    end
else
    if isempty(mask)
        % Register image to atlas with or without corresponding points
        if ~use_points
            [reg_params.img_to_atlas,reg_img]=elastix(mov_img,atlas_img,...
                outputDir,parameter_paths,'threads',[]);
        else
            [reg_params.img_to_atlas,reg_img]=elastix(mov_img,atlas_img,...
                outputDir,parameter_paths,'fp',atlas_points,'mp',mov_points,'threads',[]);
        end
    else
        % Register image to atlas with or without corresponding points +
        % with a mask
        if ~use_points
            [reg_params.img_to_atlas,reg_img]=elastix(mov_img,atlas_img,...
                outputDir,parameter_paths,'mMask',mov_mask,'fMask',atlas_mask,'threads',[]);
        else
            [reg_params.img_to_atlas,reg_img]=elastix(mov_img,atlas_img,...
                outputDir,parameter_paths,'mMask',mov_mask,'fMask',atlas_mask,...
                'fp',atlas_points,'mp',mov_points,'threads',[]);
        end
    end
end

% Image sizes
size_mov = size(mov_img);
size_atlas = size(atlas_img);

% Calculate the inverse
if isequal(calc_inverse,"true")
    fprintf('%s\t Getting inverse transformation parameters\n',datetime('now'));
    if isequal(direction, 'atlas_to_img')
        inv_direction = 'img_to_atlas';
        reg_params.atlas_to_img = get_inverse_transform_from_atlas(reg_params.img_to_atlas,...
            config, inv_direction, mov_img);
        reg_params.atlas_to_img.TransformParameters{1}.Size = size_atlas([2 1 3]);
    else
        inv_direction = 'atlas_to_img';
        reg_params.atlas_to_img = get_inverse_transform_from_atlas(reg_params.img_to_atlas,...
            config, inv_direction);
        reg_params.atlas_to_img.TransformParameters{1}.Size = size_mov([2 1 3]);
    end
end

% Remove transformed images from structure
if ~isempty(reg_params.atlas_to_img.transformedImages)
    reg_params.atlas_to_img.transformedImages = [];
end
if ~isempty(reg_params.img_to_atlas.transformedImages)
    reg_params.img_to_atlas.transformedImages = [];
end

% Remove temporary directory
rmdir(outputDir,'s')

% Check results
%imshowpair(reg_img(:,:,120),atlas_img(:,:,120))

% Save transform parameters to variables folder
reg_params.perm = perm;
reg_params.res_adj = res_adj;

% Save sizes in the registration file
reg_params.native_img_size = native_img_size([2 1 3]);
reg_params.mov_img_size = size_mov([2 1 3]);
reg_params.atlas_img_size = size_atlas([2 1 3]);

% Save registration parameters
save(fullfile(output_directory,'variables','reg_params.mat'),'reg_params')

% Write registered image
if isequal(save_registered_image, 'true')
    fprintf('%s\t Saving registered image\n',datetime('now')); 
    reg_img = uint16(reg_img);
    reg_img = imadjustn(reg_img);
    reg_img = im2uint8(reg_img);

    save_path = sprintf('%s_registered.tif',mov_img_path(1:end-4));
    imwrite(reg_img(:,:,1), save_path)
    for i = 2:size(reg_img,3)
       imwrite(reg_img(:,:,i),save_path,'WriteMode','append'); 
    end
end

end


function [mov_points,atlas_points] = load_points_from_bdv(output_directory,points_file)
% Function to load points from FIJI's Big Data Viewer
pts_path = fullfile(output_directory,'variables',points_file);
pts = readmatrix(pts_path);
pts = pts(:,3:end);

%moving x,y,z then atlas x,y,z
mov_points = pts(:,1:3);
atlas_points = pts(:,4:6);

end
