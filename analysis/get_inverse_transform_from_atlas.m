function [reg_params_inv,reg_img] = get_inverse_transform_from_atlas(reg_params, config)
% This function calculates the inverse of img_to_atlas transforms to go 
% atlas_to_img using elastix's DisplacementMagnitudePenalty metric. This
% metric doesn't give an exact inverse but an approximation, which is still
% usually pretty good. To make the approximation even more precise, you can
% decrease the grid spacing. 

downsample_factor = 0.4;

atlas_path = fullfile(config.home_path,'atlas',config.atlas_file);
atlas_img = niftiread(atlas_path);
atlas_img = imresize3(atlas_img,downsample_factor);

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

parameter_path{1} = fullfile(config.home_path,'elastix_parameter_files',...
                'atlas_registration','ElastixParameterPointsInverse.txt');

[reg_params_inv,reg_img]=elastix(atlas_img,atlas_img,outputDir,parameter_path,...
    't0',fname1,'threads',[]);

% Remove initial transform specification and reset size, spacing to match
% the moving image
reg_params_inv.TransformParameters{1}.InitialTransformParametersFileName =...
    'NoInitialTransform';

rmdir(outputDir,'s')

end