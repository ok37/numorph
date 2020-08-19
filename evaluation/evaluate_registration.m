function evaluate_registration(index)

load 'TCa_variables.mat'
config = load(fullfile(home_path,'TCa_variables.mat'));

%% Read Filename Information
path_cell = {};
    
% Load images from resampled directory
fprintf('%s\t Reading image filename information from resampled directory \n',datetime('now'))
if exist(fullfile(output_directory,'resampled'),'dir') == 7
    path_cell{1} = dir(fullfile(output_directory,'resampled'));
    location = 'resampled';
    path_table_resampled = path_to_table(path_cell,markers,channel_num,location,sample_name);
else
    msg = sprintf('%s\t Cannot locate resampled directory in the specified image directory',datetime('now'));
    error(msg);
end

mov_img_path = path_table_resampled.file{1};
config.save_registered_image = 'false';

%% Evaluate Registration Accuracy
% This registers image to atlas and would apply to registration parameters
% to the true positive, manually traced mask. Not totally representitive of
% true DICE but faster. To get true overlap of mask on native sampple, need
% to calculate inverse transformation parameters and apply these to the
% annotation image. Then measure DICE on manually traced mask without
% transforming it. For points registration, specify how many points to test
% in num_points.

if index(1)
num_points = 200;
similarity = zeros(1,length(num_points));
for n = 1:length(num_points)
    reg_params = register_to_atlas(mov_img_path, config, num_points(n));

    % Find true positive, manually traced mask. This is at 10um resolution
    mask_location = fullfile(output_directory,'resampled');
    files = dir(mask_location);
    file_idx = arrayfun(@(s) contains(s.name,'_mask.tif'), files);
        
    % Load manually traced mask and resize
    I_true = loadtiff(fullfile(files(file_idx).folder,files(file_idx).name));
    for i = 1:length(reg_params.TransformParameters)
        reg_params.TransformParameters{i}.FinalBSplineInterpolationOrder = 0;
    end
    I_true = transformix(I_true,reg_params,[1 1 1], []);
    
    % Generate annotation mask for cortical structures
    I_mask = gen_mask(hemisphere, structures_of_interest, reg_params.res_adj);    
    I_mask = imresize3(I_mask,0.4,'Method','nearest');
    
    % Calculate DICE
    similarity(n) = dice(logical(I_true),logical(I_mask));
    disp(similarity)
    
    niftiwrite(uint16(I_mask),'25um_mask.nii');
end
end
%% Create montage
% Create montage of registration results
if index(2)
n_slices = size(I_true,3);
pos = round(linspace(1,n_slices,6));

obj = cell(1,length(pos)-2);
for i = 2:length(pos)-1
    obj{i-1} = imfuse(logical(I_true(:,:,pos(i))),logical(I_mask(:,:,pos(i))));    
end
    
m = montage(obj,'Size',[1,4]);  
end
