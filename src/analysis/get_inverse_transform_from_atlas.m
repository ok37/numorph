function [reg_params_inv,reg_img] = get_inverse_transform_from_atlas(config, reg_params, direction, mov_img)
% This function calculates the inverse of img_to_atlas transforms to go 
% atlas_to_img using elastix's DisplacementMagnitudePenalty metric. This
% metric doesn't give an exact inverse but an approximation, which is still
% usually pretty good. To make the approximation even more precise, you can
% decrease the grid spacing. 

downsample_factor = 0.4;

% Note: this is the default parameter file for calculating the inverse
parameter_path{1} = fullfile(config.home_path,'elastix_parameter_files',...
                'atlas_registration','ElastixParameterPointsInverse.txt');

% Load registration parameters if not provided
if nargin<2 || isempty(reg_params)
    reg_params = load(fullfile(config.output_directory,'variables','reg_params.mat'),'reg_params');
    reg_params = reg_params.img_to_atlas;
end

if nargin<3
    error("Specify direction as atlas_to_img or img_to_atlas")
end

if isequal(direction,'img_to_atlas') && nargin<3
    error("Reference image required to register from img_to_atlas")
end

% Load ARA if going atlas_to_img
if isequal(direction,'atlas_to_img')
    if isfile(config.atlas_file)
        atlas_path = config.atlas_file;
    else
        atlas_path = fullfile(config.home_path,'supplementary_data',config.atlas_file); 
        if ~isfile(atlas_path)
            error("Could not locate Allen Reference Atlas .nii file specified")
        end
    end
    mov_img = niftiread(atlas_path);
    mov_img = imresize3(mov_img,downsample_factor);
end

% Create temporary directory for saving images
outputDir = fullfile(config.home_path,'elastix_parameter_files',...
    sprintf('%d_%d',yyyymmdd(datetime),randi(1E6,1)));
if ~exist(outputDir,'dir')
    mkdir(outputDir)
end

% To calculate inverse, only the last transform (B-spline) is needed. But
% any initial rigid transforms need to be saved as text files and specified
% in the parameters structure
if length(reg_params.TransformParameters) > 1
    for i = 1:length(reg_params.TransformParameters)-1
        reg_params.TransformParameters{i}.InitialTransformParametersFileName =...
            'NoInitialTransform';
        fname = fullfile(outputDir,sprintf('init_tform%d.txt',i));
        elastix_paramStruct2txt(fname,reg_params.TransformParameters{i})
        reg_params.TransformParameters{i+1}.InitialTransformParametersFileName = fname;
    end
end

fname1 = fullfile(outputDir,sprintf('tform_%s.txt','final'));
elastix_paramStruct2txt(fname1,reg_params.TransformParameters{end})

[reg_params_inv,reg_img]=elastix(mov_img,mov_img,outputDir,parameter_path,...
    't0',fname1,'threads',[]);

% Remove initial transform specification
reg_params_inv.TransformParameters{1}.InitialTransformParametersFileName =...
    'NoInitialTransform';

rmdir(outputDir,'s')

end