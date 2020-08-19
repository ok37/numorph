function [t_adj, lowerThresh, upperThresh, y_adj, flatfield, darkfield] = measure_images_slice(stack, config, channel_idx)
%--------------------------------------------------------------------------
% Calculate various intensity adjustments including tile differences,
% light-sheet width correction, flatfield + darkfield shading correction.
% Old version - uses slices instead of stacks
%--------------------------------------------------------------------------
low_prct = 10;   % Low percentile for sampling background pixels
high_prct = 90; % High percentile for sampling bright pixels
pads = 0.20;     % Crop this fraction of image from the sides of images
image_sampling = 0.05;   % Fraction of all images to sample

% Load config parameters
overlap = config.overlap;
ls_width = config.ls_width(channel_idx);
single_sheet = config.single_sheet;
resolution = config.resolution;
laser_y_displacement = config.laser_y_displacement(channel_idx);

% Count number of images and measure image dimensions
x_tiles = length(unique(stack.x));
y_tiles = length(unique(stack.y));
nb_images = height(stack)/(x_tiles*y_tiles);

% Get image positions
s = round(image_sampling*nb_images);
img_range = round(linspace(1,nb_images,s));

% Read image size and get padding
tempI = imread(stack.file{1});
[nrows, ncols] = size(tempI);
pad_h = round(pads*ncols*overlap);
pad_v = round(pads*nrows*overlap);

%Determine overlap region
%Horizontal image overlap regions
overlap_min_c = 1:round(ncols*overlap)-pad_h;
overlap_max_c = ncols-round(ncols*overlap)+1:ncols-pad_h;
overlap_min_h = {[1,nrows],[1,overlap_min_c(end)]};
overlap_max_h = {[1,nrows],[overlap_max_c(1),overlap_max_c(end)]};

% Vertical image overlap regions
overlap_min_r = 1:round(nrows*overlap)-pad_v;
overlap_max_r = nrows-round(nrows*overlap)+1:nrows-pad_v;
overlap_min_v = {[1,overlap_min_r(end)],[1,ncols]};
overlap_max_v = {[overlap_max_r(1),overlap_max_r(end)],[1,ncols]};

% Calculate and adjust for laser width
if isequal(config.adjust_ls_width,"true")
    fprintf('%s\t Adjusting For Laser Width \n',datetime('now'));    
    y_adj = adjust_intensity_measured(single_sheet,tempI,ls_width,resolution,laser_y_displacement)';
else
    fprintf('%s\t No Laser Width Adjustment \n',datetime('now'));    
    y_adj = ones(nrows,1);
end  

% Calculate shading correction using BaSiC
if isequal(config.shading_correction,"true")
    [flatfield, darkfield]  = estimate_flatfield(stack, config);
    flatfield = single(flatfield);
    darkfield = single(darkfield);
else
    flatfield = ones(nrows,ncols,'single'); darkfield = ones(nrows,ncols,'single');
end

% Measure pairwise horizontal intensity differences
fprintf('%s\t Measuring Between Tile Differences Horizontally \n',datetime('now'));    
h_matrix1 = ones(y_tiles,x_tiles);
h_matrix2 = ones(y_tiles,x_tiles);

% Initialize measurement vectors
p_low = zeros(1,y_tiles*(x_tiles-1)+(y_tiles-1)*x_tiles);
p_high = zeros(1,y_tiles*(x_tiles-1)+(y_tiles-1)*x_tiles);
stdev = zeros(1,y_tiles*(x_tiles-1)+(y_tiles-1)*x_tiles);

%Measure pairwise horizontal intensity differences
h_matrix = ones(y_tiles,x_tiles);
img_idx = 1;
for i = 1:y_tiles
    for j = 1:x_tiles-1
        mean_diff = zeros(1,length(img_range));
        l = 1;
        for k = img_range
            %Read image regions where tiles should overlap
            file_left = stack(stack.y == i & stack.x == j & stack.z == k,:);
            file_right = stack(stack.y == i & stack.x == j+1 & stack.z == k,:);

            I_left = double(imread(file_left.file{1},'PixelRegion',overlap_max_h))*h_matrix(i,j).*y_adj;
            I_right = double(imread(file_right.file{1},'PixelRegion',overlap_min_h))*h_matrix(i,j+1).*y_adj;     

            %Difference in average intensity indicates between tile differences
            mean_diff(l) = prctile(I_left(:),high_prct)/prctile(I_right(:),high_prct); 

            %Measure 1 percentile of all pixels. This gives rough background
            %and upper intensity
            p_low(img_idx) = mean([prctile(I_left(:),1) prctile(I_right(:),1)]);
            p_high(img_idx) = max(max([I_left(:) I_right(:)]));
            stdev(img_idx) = std2([I_left I_right]);

            %Update vector indices
            img_idx = img_idx+1;
            l = l+1;
        end
        h_matrix(i,j+1) = median(mean_diff);
    end
end
h_matrix = h_matrix/mean2(h_matrix);

%Measure pairwise vertical intensity differences
fprintf(strcat(char(datetime('now')),'\t Measuring Between Tile Differences Vertically\n'));
v_matrix = ones(y_tiles,1);
h_mean = mean(h_matrix,2);

for i = 1:y_tiles-1
    mean_diff = zeros(1,length(img_range));
    l = 1;
    for k = img_range
        %Read image regions where tiles should overlap
        file_top = stack(stack.y == i & stack.z == k,:);
        file_bottom = stack(stack.y == i+1 & stack.z == k,:);
        
        for kk = 1:height(file_top)
        
        I_top = double(imread(file_top.file{kk},'PixelRegion',overlap_max_v))*v_matrix(i);
        I_bottom = double(imread(file_bottom.file{kk},'PixelRegion',overlap_min_v))*v_matrix(i+1);
        
        %Difference in average intensity indicates between tile differences
        mean_diff(l) = prctile(I_top(:),high_prct)/prctile(I_bottom(:),high_prct);
                  
        %Measure 1 percentile of all pixels. This gives rough background
        p_low(img_idx) = mean([prctile(I_top(:),1) prctile(I_bottom(:),1)]);
        p_high(img_idx) = max(max([I_top(:), I_bottom(:)]));
        stdev(img_idx) = std2([I_top I_bottom]);
            
        %Update vector indices
        img_idx = img_idx+1;
        l = l+1;
        end
    end
    v_matrix(i+1) = median(mean_diff);
end

v_matrix = v_matrix/mean2(v_matrix);

% Combine matrices by column-wise multiplicaiton of each column
adj_matrix1 = zeros(size(h_matrix));
for i = 1:size(h_matrix,1)
    adj_matrix1 = adj_matrix1 + v_matrix.*h_matrix(i,:)/size(h_matrix,1); 
end

adj_matrix2 = zeros(size(v_matrix));
for i = 1:size(v_matrix,2)
    adj_matrix2 = adj_matrix2 + v_matrix(:,i).*h_matrix/size(v_matrix,2);
end

adj_matrix = (adj_matrix1+adj_matrix2)/2;

%Calculate lower threshold
lowerThresh = round(median(p_low)+1*median(stdev));
upperThresh = min(prctile(p_high,90),65535);

% Tile adjustment
% Save tile adjustment
t_adj(:,:,1) = adj_matrix;
t_adj(:,:,2) = adj_matrix;
disp(adj_matrix)

end