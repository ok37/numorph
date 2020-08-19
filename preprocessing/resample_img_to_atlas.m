function resample_img_to_atlas(path_table, output_directory, res, resample_res)

%Create directory for resampled images
if exist(char(strcat(output_directory,'/resampled')),'dir') ~= 7
  mkdir(char(strcat(output_directory,'/resampled')))
end

nb_images = height(path_table);

%Measure image dimensions for high resolution image
tempI = imread(path_table.file{1});
[img_height,img_width] = size(tempI);

%Calculate image dimensions for target resolution
re_v = resample_res./res;
re_height = round(img_height/re_v(1));
re_width = round(img_width/re_v(2));
re_z = round(nb_images/re_v(3));

%Intialize matrix for resampled image
re_I = zeros(re_height,re_width,re_z,'uint16');

% Resample first only in x,y to decrease memory demands
fprintf(strcat(char(datetime('now')),'\t Interpolating in xy\n'))
for i = 1:nb_images
    I = loadtiff(path_table.file{i});
    re_I(:,:,i) = imresize(I,[re_height,re_width]);
end

% Then resample in z 
fprintf(strcat(char(datetime('now')),'\t Interpolating in z\n'))
re_I = imresize3(re_I,[re_height,re_width,re_z]);
re_I = uint16(re_I);

% Save as nii series
save_directory = strcat(output_directory,'/resampled/');
filename = strcat(path_table.sample_name{1},'_',path_table.channel_num(1),'_',path_table.markers(1),'_resampled.nii');
niftiwrite(re_I, strcat(save_directory,filename))

end