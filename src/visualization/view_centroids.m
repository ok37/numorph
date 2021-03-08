function I = view_centroids(config, index, path_table, centroids_file, adjust)
%--------------------------------------------------------------------------
% View centroids overlayed onto a selected z position.
%--------------------------------------------------------------------------
% Usage:  
% I = view_centroids(config, index, path_table, centroids_file)
% 
%--------------------------------------------------------------------------
% Inputs:
% config: Configuration structure for NM_analyze.
%
% index: (1x2 int) Indexes specifying image channel (1) and z position (2).
%
% path_table: Path table specifying image filenames. (Optional: will
% attempt to reload from config if not specified)
%
% centroids_file: (string) Filename specifying specific centroids csv file.
% (Default: load most recently generated centroids file) 
%
% adjust: (logical) Rescale image intensity. (Default: true)
%
%--------------------------------------------------------------------------
% Output:
% I: (3 channel uint8 image) Output image with centroids overlayed.
%
%--------------------------------------------------------------------------

if length(index) ~= 2
    error("Specify channel and z position index as 1x2 vector")
end

if nargin<3
    % Load path table based on config as in NM_analyze
    % Generate table containing image information
    path_table = path_to_table(config);

    % If markers ignored, add these from raw
    if ~all(ismember(config.markers,unique(path_table.markers)))
        idx = find(~ismember(config.markers,unique(path_table.markers)));
        for i = 1:length(idx)
            config2 = config;
            config2.markers = config.markers(idx(i));
            config2.channel_num = config.channel_num(idx(i));
            try 
                path_table = vertcat(path_table,path_to_table(config2,'raw',false));
                path_table.channel_num = arrayfun(@(s) find(s == config.markers),path_table.markers);
            catch
                warning("Could not load marker %s ignored from processing.",config.markers(idx(i)));
            end
        end
        clear config2
    end
end

if nargin<4
    % Find most recent centroids file
    files = dir(config.output_directory);
    files = files(arrayfun(@(s) contains(s.name,'_centroids') &&...
        endsWith(s.name,'.csv'), files));

    if length(files)>1
        [~,idx] = sort([files.datenum],'descend');
        files = files(idx);
    end
    centroids = readmatrix(fullfile(config.output_directory,files(1).name));
else
    cen_path = fullfile(config.output_directory,centroids_file);
    if isfile(cenpath)
        centroids = readmatrix(cen_path);
    else
        error("Could not locate specified centroids file")
    end
end

if nargin < 5
    adjust = true;
end

% Read image from path_table
I = read_img(path_table,index);

if islogical(adjust) && adjust
    I = imadjust(I);
elseif isnumeric(adjust)
    I = imadjust(I,[adjust(1) adjust(2)]);
end

% Subset centroids
z_position = index(2)-1:index(2)+1;
idx = any(centroids(:,3) == z_position,2);
centroids = centroids(idx,:);

% Create binary mask
BW = zeros(size(I));
cen_mid = centroids(centroids(:,3) == z_position(2),:);
cen_mid = sub2ind(size(BW),cen_mid(:,1),cen_mid(:,2));
BW(cen_mid) = 1;

% Make a cross for cells in current slice
se = strel([0 1 0;1 1 1;0 1 0]);
BW = imdilate(BW,se);

% Add points to centroids above and below
cen_up = centroids(centroids(:,3) == z_position(1),:);
cen_up = sub2ind(size(BW),cen_up(:,1),cen_up(:,2));
BW(cen_up) = 1;
cen_down = centroids(centroids(:,3) == z_position(3),:);
cen_down = sub2ind(size(BW),cen_down(:,1),cen_down(:,2));
BW(cen_down) = 1;

% Overlay mask
I = imoverlay(I,BW,'r');

if nargout<1
    imshow(I)
end

end