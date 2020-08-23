function resample_img(config, resample_res, markers, file_type)
%--------------------------------------------------------------------------
% Resample mage to specified resolution. Only requires NM_analysis config 
% structure, desired resolution, and image filetype.
%--------------------------------------------------------------------------

if length(resample_res) == 1
    resample_res = repmat(resample_res,1,3);
end

if nargin<3 || isempty(markers)
    markers = config.markers;
end

if nargin<4
    file_type = '.nii';
end

% Create directory for resampled images
if exist(fullfile(config.output_directory,'resampled'),'dir') ~= 7
  mkdir(fullfile(config.output_directory,'resampled'))
end

% Create load table with file info
path_cell{1} = dir(config.img_directory);
if isequal(fullfile(config.output_directory,'stitched'),config.img_directory)
    fprintf('%s\t Reading image filename information from stitched directory \n',datetime('now'))
    location = 'stitched';
    path_table = path_to_table(path_cell,location,[],[],config.sample_name);
end

for i = 1:length(markers)
    % Subset marker
    path_sub = path_table(path_table.markers == markers(i),:);
    path_sub = path_sub(1:10,:);
    nb_images = height(path_sub);
    
    %Measure image dimensions for high resolution image
    tempI = imread(path_sub.file{1});
    [img_height,img_width] = size(tempI);

    %Calculate image dimensions for target resolution
    re_v = resample_res./config.resolution;
    re_height = round(img_height/re_v(1));
    re_width = round(img_width/re_v(2));
    re_z = round(nb_images/re_v(3));

    %Intialize matrix for resampled image
    re_I = zeros(re_height,re_width,re_z,'uint16');

    % Resample first only in x,y to decrease memory demands
    fprintf('%s\t Resampling in xy \n',datetime('now'))
    for j = 1:nb_images
        I = loadtiff(path_sub.file{j});
        re_I(:,:,j) = imresize(I,[re_height,re_width]);
    end

    % Then resample in z 
    fprintf('%s\t Resampling in z \n',datetime('now'))
    re_I = imresize3(re_I,[re_height,re_width,re_z]);
    re_I = uint16(re_I);

    % Save as nii series
    save_directory = fullfile(config.output_directory,'resampled');
    channel_idx = find(config.markers == markers(i));
    
    if isequal(file_type,'nii')
        filename = sprintf("%s_C%d_%s_resampled_%d.nii",config.sample_name,channel_idx,config.markers(channel_idx),resample_res(1));
        niftiwrite(re_I, fullfile(save_directory,filename))
    elseif isequal(file_type,'.tif') || isequal(file_type,'.tiff')
        filename = sprintf("%s_C%d_%s_resampled_%d.tif",config.sample_name,channel_idx,config.markers(channel_idx),resample_res(1));
        saveastiff(re_I, char(fullfile(save_directory,filename)))
    end
end

end