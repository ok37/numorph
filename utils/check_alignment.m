function [img] = check_alignment(config, ranges, markers, spacing)
%--------------------------------------------------------------------------
% Check image alignment of channels
%--------------------------------------------------------------------------
% Inputs:
% config - config structure from NM_process.
%
% ranges - 1x2 or 1x3 cell array containing tile row and column. 3d value
% can be a range of slices in the slices (e.g. 1:100).
%
% markers - (optional) string array containing which markers to produce
% images for. Can also also specify as numeric for marker number.
%
% spacing - (optional) 1x3 double for the amount of downsampling for each
% dimension.

if ~iscell(ranges)
    error("Second input should be 1x2 cell array containing tile row, column position. "+...
        "Optional 3rd cell array input to subset a range of tile slices. \n")
end

if nargin<3
    markers = config.markers;
elseif isnumeric(markers)
    markers = config.markers(markers);
end

if nargin<4 || isempty(spacing)
    spacing = [5,5,10];
end

fprintf("Loading image information \n")

% Load images from aligned directory
img_directory = fullfile(config.output_directory,'aligned');
path_cell{1} = dir(img_directory);
location = "aligned";
path_table = path_to_table(path_cell,location,[],[]);

% Subset x,y positions
path_table = path_table(ismember(path_table.markers,markers),:);
path_table = path_table(path_table.y == ranges{1} & ...
    path_table.x == ranges{2},:);
total_z = height(path_table)/(length(markers));

if length(ranges) == 2
    z_min = 1;
    z_max = total_z;
elseif length(ranges{3}) == 1
    z_min = 1;
    z_max = max(ranges{3});
else
    z_min = min(ranges{3});
    z_max = max(ranges{3});
end

z_pos = z_min:spacing(3):z_max;
path_table = path_table(ismember(path_table.z,z_pos),:);

tempI = imread(path_table.file{1});
[nrows,ncols] = size(tempI);
nrows_d = round(nrows/spacing(1));
ncols_d = round(ncols/spacing(2));

img = zeros(nrows_d,ncols_d,length(z_pos),length(markers),'uint16');
for i = 1:length(markers)
    fprintf("Reading marker %s \n",markers(i)) 
    path_sub = path_table(path_table.markers == markers(i),:);
    for j = 1:length(z_pos)
       I = imread(path_sub.file{j});
       img(:,:,j,i) = imresize(I,[nrows_d,ncols_d]);
    end
    fname = fullfile(config.output_directory,'samples',sprintf('%s_%d_%d.tif',...
        markers(i),ranges{1},ranges{2}));
    options.overwrite = true;
    saveastiff(squeeze(img(:,:,:,i)), char(fname),options)
end

end