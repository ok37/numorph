function [reg_params_inv,reg_img] = get_inverse_registration_transform(reg_params,...
    mov_img_path, home_path)
% This function calculate the inverse transform using elastix's
% DisplacementMagnitudePenalty metric. Everything in the inverse transform
% parameter text file should match the original transform except for this
% metric

% Read moving image
mov_img = niftiread(mov_img_path);
if isequal(class(mov_img),'double')
    mov_img = flip(mov_img,1);
    mov_img = imrotate(mov_img,-90);
else
   mov_img = double(mov_img); 
end

% Resize image
mov_img = imresize3(mov_img,0.4);

% Adjust intensity to match those in register_to_atlas.m
upperThresh = stretchlim(mov_img(:));
mov_img = double(imadjustn(uint16(mov_img),[0 upperThresh(2)],[],0.9));

% Create temporary directory for saving images
i = 33;
outputDir = fullfile(home_path,'elastix_parameter_files',sprintf('tmp%d',i));
if ~exist(outputDir,'dir')
    mkdir(outputDir)
end

% Location of inverse transform file. This should match 
parameter_path{1} = fullfile(home_path,'elastix_parameter_files',...
                'atlas_registration','ElastixParameterPointsInverse.txt');

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

fname1 = fullfile(outputDir,sprintf('hello_tform_%s.txt','final'));
elastix_paramStruct2txt(fname1,reg_params.TransformParameters{end})

% Run the registration using the image of interest as the fixed and moving
% image
[reg_params_inv,reg_img]=elastix(mov_img,mov_img,outputDir,parameter_path,...
    't0',fname1,'threads',[]);

%imshowpair(uint16(mov_img(:,:,100)),uint16(reg_img(:,:,100)))

% Delete temporary directory
%rmdir(outputDir,'s')

% Remove initial transform specification and reset size, spacing to match
% the moving image
reg_params_inv.TransformParameters{1}.InitialTransformParametersFileName =...
    'NoInitialTransform';
%img_size = size(mov_img);
%img_size = img_size([2 1 3]);
%reg_params_inv.TransformParameters{1}.Spacing = img_size ./...
%    reg_params.TransformParameters{end}.Size;

end