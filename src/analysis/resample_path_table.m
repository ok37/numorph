function resampled_image = resample_path_table(path_table, config, resolution, resample_res, channels)
%--------------------------------------------------------------------------
% Resample image to specific resolution. 
%--------------------------------------------------------------------------
% Inputs:
% path_table - table of image paths. Images should contain only 1 tile.
%
% config - config structure from NM_analyze. Optional: can leave empty and
% specify resoltion and resample_res.
%
% resolution - [x,y,z] specifying input image resolution. Can also send as
% cell array of resolutions if different for each channel. 
%
% resample_res - [x,y,z] resolution to resample to.
%
% channels - (default: length) index specifying which channels to resample.
%--------------------------------------------------------------------------

% Must specify resample resolution somewhere
if isempty(config) && nargin<3
    error("Must provide config structure or specify resolutions.")
end

% Default resample resolution
if nargin<4
    resample_res = config.resample_resolution;
end

% Default do all channels in path_table
if nargin<5
    channels = unique(path_table.channel_num);
end
    
% Check resolution
if nargin>3 && ~iscell(resolution)
    resolution = repmat({resolution},1,length(channels));
else
    resolution = config.resolution;
end

% Create resampled directory
resample_path = fullfile(config.output_directory,'resampled');
if exist(resample_path,'dir') ~= 7
    mkdir(resample_path)
end

% Perform resampling for 1 or multiple channels
nchannels = length(channels);
resampled_image = cell(1,nchannels);
for i = 1:nchannels
    path_sub = path_table(path_table.channel_num == channels(i),:);
    nb_images = height(path_sub);

    % Measure image dimensions for high resolution image
    tempI = imread(path_sub.file{1});
    [nrows,ncols] = size(tempI);

    % Calculate image dimensions for target resolution
    re_v = resample_res./resolution{i};
    re_height = round(nrows/re_v(1));
    re_width = round(ncols/re_v(2));
    re_z = round(nb_images/re_v(3));

    % Intialize matrix for resampled image
    re_I = zeros(re_height,re_width,re_z,'uint16');

    % Resample first only in x,y to decrease memory demands
    fprintf('Reading and resampling XY\n')            
    for j = 1:nb_images
        I = loadtiff(path_table.file{i});
        re_I(:,:,j) = imresize(I,[re_height,re_width]);
    end

    % Then resample in z 
    fprintf('Resampling Z\n')            
    re_I = imresize3(re_I,[re_height,re_width,re_z]);

    % Save as nii series
    marker = path_table(path_table.channel_num == channels(i),:).markers{1};
    if ~isempty(config) && isfield(config,'sample_name')
        filepath = fullfile(resample_path,sprintf('%s_C%d_%s_%d_%d_%d.nii',config.sample_name,...
            channels(i),marker,resample_res(1),resample_res(2),resample_res(3)));
    else
        filepath = fullfile(resample_path,sprintf('C%d_%s_%d_%d_%d.nii',...
            channels(i),marker,resample_res(1),resample_res(2),resample_res(3)));
    end
    niftiwrite(re_I, filepath)
    
    % Check if providing resampled image
    if nargout == 1
        resampled_image{i} = re_I;
    end
end

end