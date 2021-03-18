function reg_params = register_to_atlas(config, mov_img_path, num_points, only_inverse)
% Register images to the reference atlas using elastix via melastix
% wrapper. Note: the final registration parameters are stored in the
% reg_params variable. While it is simpler to just register the atlas to
% the image, you may, in some cases, get better accuracy by registering the
% image to the atlas and calculating the inverse parameters. 

% Defaults
high_prct = 99;                 % Max intensity adjustment percentile for images
Gamma = 0.9;                    % Gamma to apply for intensity adjustment
inverse_params = "inverse";     % Location of elastix parameter for calculating the inverse

% Unpack config variables
params = config.registration_parameters;
resample_res = config.resample_resolution;
hemisphere = config.hemisphere;
orientation = config.orientation;
output_directory = config.output_directory;
use_mask = config.mask_cerebellum_olfactory;
direction = config.direction;
calc_inverse = config.calculate_inverse;
points_file = config.points_file;
atlas_file = config.atlas_file;
perm = [0,0];
home_path = fileparts(which('NM_config.m'));

% Unless testing, use all points
if nargin == 2
    num_points = 'all';
end

if nargin<4
    only_inverse = false;
end

% Location of parameters
if isempty(params)
    if isempty(points_file)
        params = "default";
    else
        params = "points";
    end
end

% Type of atlas
if contains(atlas_file,'nissl')
    atlas_type = 'nissl';
elseif contains(atlas_file,'average')
    atlas_type = 'average';
else
    atlas_type = 'nissl';
end
    
% Convert to cell
if ~iscell(mov_img_path)
    mov_img_path = {mov_img_path};
end

% Check atlas files
if iscell(atlas_file)
    atlas_path = cellfun(@(s) fullfile(config.home_path,'data','atlas',s),atlas_file,'UniformOutput',false);
else
    atlas_path = arrayfun(@(s) {char(fullfile(config.home_path,'data','atlas',s))},atlas_file);
end
assert(all(isfile(atlas_path)), "Could not locate Allen Reference Atlas .nii file specified")

% Load atlas + moving image 
fprintf('Reading resampled and atlas images\n'); 

% MATLAB will usually rotate and flip axis dimensions based on the data
% type in the niftifile
mov_img_array = cell(1,length(mov_img_path));
for i = 1:length(mov_img_path)
    mov_img = niftiread(mov_img_path{i});
    if isequal(class(mov_img),'double')
        mov_img = imrotate(mov_img,90);
        mov_img = flip(mov_img,1);
    else
       mov_img = double(mov_img); 
    end
    mov_img_array{i} = mov_img;
end

% Load ARA
atlas_img_array = cell(1,length(atlas_path));
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
            warning('Using lateral orientation on whole brain image\n');
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
    atlas_img_array{i} = double(atlas_img);
end

% Resize
native_img_size = size(mov_img_array{1});
atlas_res = cellfun(@(s) repmat(str2double(regexp(s,'\d*','Match')),1,3),atlas_file,'UniformOutput',false);
% Rule: all atlas file must be at the same resolution
assert(all(atlas_res{1} == atlas_res{end}),"All loaded atlas file must be at the same resolution")
res_adj = atlas_res{1}./resample_res;
if res_adj(1) > 1
    % Atlas is at a lower resolution. Downsample image
    mov_img_array = cellfun(@(s) imresize3(s,res_adj.*size(s)),mov_img_array);
elseif res_adj(1) < 1
    % Image is at a lower resolution. Downsample atlas
    atlas_img_array = cellfun(@(s) imresize3(s,res_adj.*size(s)),atlas_img_array);
end

% Rescale intensities
mov_img_array = cellfun(@(s) double(imadjustn(uint16(s),...
    [0,prctile(s(:),high_prct)/65535],[],Gamma)), mov_img_array,...
    'UniformOutput', false);
    
% Chain registration parameters. Parameter files are found in
% ./elastix_parameter_files/atlas_registration. Update specific parameters
% by editting these files.

if ~only_inverse
    % Paths to elastix parameter files
    param_path = fullfile(home_path,'data','elastix_parameter_files','atlas_registration',params);
    param_path = dir(param_path);
    parameter_paths = cell(1,length(params));
    if ~isempty(param_path)
        % Detect which transform is in the parameter file
        parameterSub = param_path(arrayfun(@(s) endsWith(s.name,'.txt'),param_path));
        for j = 1:length(parameterSub)
            file_path = fullfile(parameterSub(1).folder,parameterSub(j).name);
            text = textread(file_path,'%s','delimiter','\n');
            n = find(cellfun(@(s) contains(s,'(Transform '),text));
            if contains(text(n),'Translation') || contains(text(n),'Affine') || contains(text(n),'Euler')
                fprintf('\t Performing rigid registration\n');
                parameter_paths{1} = file_path;
            elseif contains(text(n),'BSpline')
                fprintf('\t Performing b-spline registration\n');
                parameter_paths{2} = file_path;
            end
        end
    else
        error("Could not locate elastix parameter folder %s.",param_path)
    end

    % Load registration points
    points = [];
    if ~isempty(points_file)
        fprintf('\t Loading points to guide registration\n');

        % Load points from BigWarpJ .csv file
        [mov_points,atlas_points] = load_points_from_bdv(output_directory, points_file,...
            atlas_type, mov_img);

        % Trim points if not using all
        if isequal(num_points,'all')
            num_points = size(mov_points,1);
        end
        
        points.mov_points = mov_points(1:num_points,:);
        points.atlas_points = atlas_points(1:num_points,:);
    end

    % Add mask if coordinates are specified
    mask = [];
    %if ~isempty(mask_coordinates)
    %    mask = zeros(size(mov_img));
    %    mask(mask_coordinates(1):mask_coordinates(2),:,:) = 1;
    %end
    if isequal(use_mask,"true")
        load(fullfile(home_path,'data','annotation_data','olf_cer.mat'),'bw_mask')
        mask = single(~bw_mask);
    end

    % Perform registration
    % Single channel, pairwise registration
    [reg_params,reg_img] = elastix_registration(atlas_img_array,mov_img_array,...
        parameter_paths,points,mask,config.output_directory,direction,config.register_channels);
else
    % Load registration parameters from file
    if isfile(reg_file) 
        load(reg_file,'reg_params')
    end
end

% Image sizes
size_mov = size(mov_img);
size_atlas = size(atlas_img);

% Calculate the inverse
if isequal(calc_inverse,"true")
    fprintf('\t Getting inverse transformation parameters\n');
    if isequal(direction, "atlas_to_image")
        inv_direction = "image_to_atlas";
        reg_params.atlas_to_image = get_inverse_transform_from_atlas(config,...
            mov_img_array{1},reg_params,inv_direction,inverse_params);
        reg_params.atlas_to_image.TransformParameters{1}.Size = size_atlas([2 1 3]);
    else
        inv_direction = "atlas_to_image";
        reg_params.atlas_to_image = get_inverse_transform_from_atlas(config,...
            atlas_img_array{1},reg_params,inv_direction,inverse_params);
        reg_params.atlas_to_image.TransformParameters{1}.Size = size_mov([2 1 3]);
    end
end

% Remove transformed images from structure
if ~isfield(reg_params.atlas_to_image,'transformedImages')
    reg_params.atlas_to_image.transformedImages = [];
end
if ~isfield(reg_params.image_to_atlas,'transformedImages')
    reg_params.image_to_atlas.transformedImages = [];
end

% Check results
%imshowpair(reg_img(:,:,120),atlas_img(:,:,120))

% Save transform parameters to variables folder
reg_params.perm = perm;
reg_params.atlas_res = atlas_res;
reg_params.mov_res = resample_res;

% Save sizes in the registration file
reg_params.native_img_size = native_img_size([2 1 3]);
reg_params.mov_img_size = size_mov([2 1 3]);
reg_params.atlas_img_size = size_atlas([2 1 3]);

% Write registered image
%if isequal(save_registered_image, 'true')
%    fprintf('\t Saving registered image\n'); 
%    reg_img = uint16(reg_img);
%    reg_img = imadjustn(reg_img);
%    reg_img = im2uint8(reg_img);
    
    % Create temporary directory for saving images
%    outputDir = fullfile(config.output_directory,'registered');
%    if ~exist(outputDir,'dir')
%        mkdir(outputDir)
%    end
    
    % Save copy of completementary image
%    if isequal(direction,"atlas_to_image")
%        % Save registered image
%        [~,name] = fileparts(mov_img_path{1});
%        save_path = fullfile(outputDir,sprintf('%s_registered.tif',name));
%        imwrite(reg_img(:,:,1), save_path)
%        for i = 2:size(reg_img,3)
%            imwrite(reg_img(:,:,i),save_path,'WriteMode','append'); 
%        end

        % Save copy of moving image
%        save_path = fullfile(outputDir,sprintf('%s_registered.tif',name));
%        imwrite(reg_img(:,:,1), save_path)
%        for i = 2:size(reg_img,3)
%            imwrite(reg_img(:,:,i),save_path,'WriteMode','append'); 
%        end
%    else
%        
%        
%    end
%end

end


function [reg_params, reg_img, status] = elastix_registration(atlas_img,mov_img,...
    parameter_paths,points,mask,outputDir,direction,reg_channels)

% Single channel, pairwise registration

% Here mask is assumed to be only on the true atlas image
if isequal(direction,'atlas_to_image') 
end

% Check if using points
use_points = false;
if ~isempty(points)
    use_points = true;
    mov_points = points.mov_points;
    atlas_points = points.atlas_points;
end

% If multiple moving channels, apply rigid registration to pre-align
if length(mov_img)>1
    fprintf('\t Pre-aligning multiple input image channels prior to registration\n');
    
    % Create temporary directory for saving images
    outputDir = fullfile(outputDir,sprintf('tmp_reg_%d',randi(1E4)));
    if ~isfolder(outputDir)
        mkdir(outputDir)
    end
    
    home_path = fileparts(which('NM_config'));
    pre_align_path = {fullfile(home_path,'data','elastix_parameter_files','atlas_registration',...
        'pre-align','ElastixParameterAffine.txt')};
    for i = 2:length(mov_img)
        [~,mov_img{i},status]=elastix(mov_img{i},mov_img{1},outputDir,pre_align_path,'threads',[]);
    end
    
    % Remove temporary directory
    rmdir(outputDir,'s')
    
    % Check status and raise error if registration failed
    if status ~= 0
        error("Errors occured during registration")
    end
    
    if length(reg_channels) <length(mov_img)
        mov_img = mov_img(2:end);
    end
end


% Create temporary directory for saving images
outputDir = fullfile(outputDir,sprintf('tmp_reg_%d',randi(1E4)));
if ~isfolder(outputDir)
    mkdir(outputDir)
end

% Perform registration
reg_params = struct('image_to_atlas',[],'atlas_to_image',[]);
if isequal(direction,'atlas_to_image')
    fprintf('\t Registering atlas to input image\n');
    if isempty(mask)
        % Register atlas to image with or without corresponding points
        if ~use_points
            [reg_params.atlas_to_image,reg_img,status]=elastix(atlas_img,mov_img,...
                outputDir,parameter_paths,'threads',[]);
        else
            [reg_params.atlas_to_image,reg_img,status]=elastix(atlas_img,mov_img,...
                outputDir,parameter_paths,'fp',points.mov_points,'mp',...
                points.atlas_points,'threads',[]);
        end
    else
        % Register atlas to image with or without corresponding points +
        % with a mask
        if ~use_points
            [reg_params.atlas_to_image,reg_img,status]=elastix(atlas_img,mov_img,...
                outputDir,parameter_paths,'mMask',mask,'threads',[]);
        else
            [reg_params.atlas_to_image,reg_img,status]=elastix(atlas_img,mov_img,...
                outputDir,parameter_paths,'fp',mov_points,'mp',atlas_points,...
                'mMask',mask,'threads',[]);
        end
        % Remove transformed image from structure
        reg_params.atlas_to_image.transformedImages = [];
    end
else
    fprintf('\t Registering input image to atlas\n');
    if isempty(mask)
        % Register image to atlas with or without corresponding points
        if ~use_points
            [reg_params.image_to_atlas,reg_img,status]=elastix(mov_img,atlas_img,...
                outputDir,parameter_paths,'threads',[]);
        else
            [reg_params.image_to_atlas,reg_img,status]=elastix(mov_img,atlas_img,...
                outputDir,parameter_paths,'fp',points.atlas_points,...
                'mp',points.mov_points,'threads',[]);
        end
    else
        % Register image to atlas with or without corresponding points +
        % with a mask
        if ~use_points
            [reg_params.image_to_atlas,reg_img,status]=elastix(mov_img,atlas_img,...
                outputDir,parameter_paths,'fMask',mask,'threads',[]);
        else
            [reg_params.image_to_atlas,reg_img,status]=elastix(mov_img,atlas_img,...
                outputDir,parameter_paths,'fMask',mask,...
                'fp',points.atlas_points,'mp',points.mov_points,'threads',[]);
        end
    end
end

% Remove temporary directory
rmdir(outputDir,'s')

% Check status and raise error if registration failed
if status ~= 0
    error("Errors occured during registration")
end

end

function [mov_points,atlas_points] = load_points_from_bdv(output_directory,points_file,atlas_type,mov_img)
% Function to load points from FIJI's Big Data Viewer

files = dir(fullfile(output_directory, '**', '*.csv'));
files = files(arrayfun(@(s) isequal(s.name,points_file),files),:);

if isempty(files)
    files = dir(fullfile(output_directory, '*.csv'));
    if isempty(files)
        fprintf("\t Generating default 50 point landmarks file for %s atlas. "+...
            "Use this as a basis for point selection using Fiji's Big Warp package\n",atlas_type);
        generate_default_points(output_directory,atlas_type)
        % Save also moving image
        save_name = fullfile(output_directory,'landmarks_img_target.tif');
        saveastiff(uint16(mov_img),char(save_name));        
        pause(5)
        error("Exiting")
    else
        error("Could not locate points file %s",points_file)
    end
elseif length(files) >1
    error("Multiple files detected with the name %s in the output directory",points_file)
end

pts_path = fullfile(files(1).folder,files(1).name);
pts = readmatrix(pts_path);
pts = pts(:,3:end);

%moving x,y,z then atlas x,y,z
mov_points = pts(:,1:3);
atlas_points = pts(:,4:6);

end


function generate_default_points(output_directory,atlas_type)
% Copy pre-defined 50 points to output directory
home_path = fileparts(which('NM_config'));
def_points_path = fullfile(home_path,'data','elastix_parameter_files','atlas_registration','landmarks.csv');
m = readtable(def_points_path);

m2 = m(:,1:2);
if isequal(atlas_type,'nissl')
    m2{:,3:5} = m{:,6:8}; 
    m2 = horzcat(m2,array2table(repmat("Infinity",size(m,1),3),'VariableNames',{'a','b','c'}));
    
    % Move atlas file over to output directory
    atlas_file = fullfile(home_path,'data','atlas','ara_nissl_25.nii');
    img = niftiread(atlas_file);
    
    save_name = fullfile(output_directory,'landmarks_nissl_moving.tif');
    options.message = false;
    options.overwrite = true;
    saveastiff(img,char(save_name),options);
    save_name = 'landmarks_50_nissl.csv';
    
else
    m2{:,3:5} = m{:,3:5};
    m2 = horzcat(m2,array2table(repmat("Infinity",size(m,1),3),'VariableNames',{'a','b','c'}));
    
    % Move atlas file over to output directory
    atlas_file = fullfile(home_path,'data','atlas','average_template_25.nii');
    img = niftiread(atlas_file);
    
    save_name = fullfile(output_directory,'landmarks_average_moving.tif');
    options.message = false;
    options.overwrite = true;
    saveastiff(img,char(save_name),options);
    save_name = 'landmarks_50_average.csv';
end

save_path = fullfile(output_directory,save_name);
writetable(m2,save_path,'WriteVariableNames',0)

end