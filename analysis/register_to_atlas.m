function reg_params = register_to_atlas(mov_img_path, config, num_points, mov_mask, atlas_mask)
% Register images to the reference atlas using elastix via melastix
% wrapper. Note: the final registration parameters are stored in the
% reg_params variable. While it is simpler to just register the atlas to
% the image, you may, in some cases, get better accuracy by registering the
% image to the atlas and calculating the inverse parameters. This function
% will default to automatically measure the inverse parameters in this
% case, unless 'default_calc_inverse' is set to false.

default_calc_inverse = "true";

% Unpack config variables
resample_res = config.resample_res;
registration_method = config.registration_method;
home_path = config.home_path;
hemisphere = config.hemisphere;
orientation = config.orientation;
output_directory = config.output_directory;
mask_coordinates = config.mask_coordinates;
direction = config.direction;
points_file = config.points_file;
save_registered_image = config.save_registered_image;
atlas_file = config.atlas_file;
perm = [0,0];

% Unless testing, use all points
if nargin == 2
    num_points = 'all';
end

% Load atlas + moving image 
fprintf('%s\t Reading resampled and atlas images\n',datetime('now')); 

% MATLAB will usually rotate and flip axis dimensions based on the data
% type in the niftifile
mov_img = niftiread(mov_img_path);
if isequal(class(mov_img),'double')
    mov_img = flip(mov_img,1);
    mov_img = imrotate(mov_img,-90);
else
   mov_img = double(mov_img); 
end

atlas_path = fullfile(home_path,'atlas',atlas_file);
atlas_img = niftiread(atlas_path);

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

% Resize
native_img_size = size(mov_img);
atlas_res = repmat(str2double(regexp(atlas_file,'\d*','Match')),1,3);
res_adj = atlas_res./resample_res;
if res_adj(1) > 1
    % Atlas is at a lower resolution. Downsample image
    new_img_size = res_adj.*size(mov_img);
    mov_img = imresize3(mov_img,new_img_size);
elseif res_adj(1) < 1
    % Image is at a lower resolution. Downsample atlas
    new_atlas_size = res_adj.*size(atlas_img);
    atlas_img = imresize3(atlas_img,new_atlas_size);
end

% Rescale intensities
upperThresh = double(prctile(mov_img(:),99))/65535;
mov_img = double(imadjustn(uint16(mov_img),[0 upperThresh],[],0.9));

mov_img = double(imresize3(mov_img,0.4));
atlas_img = double(imresize3(atlas_img,0.4));

% Add mask if coordinates are specified
mask = [];
if ~isempty(mask_coordinates)
    mask = zeros(size(mov_img));
    mask(mask_coordinates(1):mask_coordinates(2),:,:) = 1;
end
    
% Chain registration parameters. Parameter files are found in
% ./elastix_parameter_files/atlas_registration. Update specific parameters
% by editting these files.
parameter_path = cell(1,length(registration_method));
use_points = 'false';
for i = 1:length(registration_method)
    reg_type = registration_method(i);
    switch reg_type
        case 'a'
            fprintf('%s\t Performing affine registration\n',datetime('now'));
            
            parameter_path{i} = fullfile(home_path,'elastix_parameter_files',...
                'atlas_registration','ElastixParameterAffine.txt');
        case 'b'
            fprintf('%s\t Performing b-spline registration\n',datetime('now')); 
            
            parameter_path{i} = fullfile(home_path,'elastix_parameter_files',...
                'atlas_registration','ElastixParameterBSpline.txt');

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
            use_points = 'true';                        
            
            parameter_path{1} = fullfile(home_path,'elastix_parameter_files',...
                'atlas_registration','ElastixParameterAffinePoints.txt');
            
            parameter_path{2} = fullfile(home_path,'elastix_parameter_files',...
                'atlas_registration','ElastixParameterPoints.txt');
            break
    end
end

% Create temporary directory for saving images
outputDir = fullfile(home_path,'elastix_parameter_files',...
    sprintf('%d_%d',yyyymmdd(datetime),randi(1E6,1)));
if ~exist(outputDir,'dir')
    mkdir(outputDir)
end

reg_params = struct('img_to_atlas',[],'atlas_to_img',[]);

% Perform registration
if isequal(direction,'atlas_to_img')
    if isempty(mask)
        % Register atlas to image with or without corresponding points
        if isequal(use_points,'false')
            [reg_params.atlas_to_img,reg_img]=elastix(atlas_img,mov_img,...
                outputDir,parameter_path,'threads',[]);
        else
            [reg_params.atlas_to_img,reg_img]=elastix(atlas_img,mov_img,...
                outputDir,parameter_path,'fp',mov_points,'mp',atlas_points,'threads',[]);
        end
    else
        % Register atlas to image with or without corresponding points +
        % with a mask
        if isequal(use_points,'false')
            [reg_params.atlas_to_img,reg_img]=elastix(atlas_img,mov_img,...
                outputDir,parameter_path,'mMask',mov_mask,'fMask',atlas_mask,'threads',[]);
        else
            [reg_params.atlas_to_img,reg_img]=elastix(atlas_img,mov_img,...
                outputDir,parameter_path,'fp',mov_points,'mp',atlas_points,...
                'mMask',mov_mask,'fMask',atlas_mask,'threads',[]);
        end
        % Remove transformed image from structure
        reg_params.atlas_to_img.transformedImages = [];
    end
elseif isequal(direction,'img_to_atlas')
    if isempty(mask)
        % Register image to atlas with or without corresponding points
        if isequal(use_points,'false')
            [reg_params.img_to_atlas,reg_img]=elastix(mov_img,atlas_img,...
                outputDir,parameter_path,'threads',[]);
        else
            [reg_params.img_to_atlas,reg_img]=elastix(mov_img,atlas_img,...
                outputDir,parameter_path,'fp',atlas_points,'mp',mov_points,'threads',[]);
        end
    else
        % Register image to atlas with or without corresponding points +
        % with a mask
        if isequal(use_points,'false')
            [reg_params.img_to_atlas,reg_img]=elastix(mov_img,atlas_img,...
                outputDir,parameter_path,'mMask',mov_mask,'fMask',atlas_mask,'threads',[]);
        else
            [reg_params.img_to_atlas,reg_img]=elastix(mov_img,atlas_img,...
                outputDir,parameter_path,'mMask',mov_mask,'fMask',atlas_mask,...
                'fp',atlas_points,'mp',mov_points,'threads',[]);
        end
    end
    % Calculate the inverse
    if isequal(default_calc_inverse,"true")
        fprintf('%s\t Getting inverse transformation parameters\n',datetime('now')); 
        [reg_params.atlas_to_img, reg_img] = get_inverse_transform_from_atlas(reg_params.img_to_atlas,...
            config);
        reg_params.atlas_to_img.TransformParameters{1}.Size = ...
            reg_params.img_to_atlas.TransformParameters{1}.Size;
    end
    % Remove transformed images from structure
    reg_params.atlas_to_img.transformedImages = [];
    reg_params.img_to_atlas.transformedImages = [];
else
    error('%s\t Incorrect direction specified \n',string(datetime('now')))
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
size1 = size(mov_img);
reg_params.mov_img_size = size1([2 1 3]);
size2 = size(atlas_img);
reg_params.atlas_img_size = size2([2 1 3]);

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

function save_points_to_bdv(mov_points,ref_points,output_directory)
% Function to save updated points that can be loaded into FIJI's Big Data
% Viewer

pts_filename = fullfile(output_directory,'variables','native_points_bdv.csv');

n_points = 1:size(mov_points,1);
pts_name =  arrayfun(@(s) sprintf('Pt-%d',s),n_points,'UniformOutput',false)';
pts_active = repmat('TRUE',size(mov_points,1),1);

pts_table = table(pts_name,pts_active,mov_points,ref_points);
writetable(pts_table,pts_filename,'WriteVariableNames',0)
end
